function [qTrue, omegaTrue, torques] = simulateAttitude(t, rECI, vECI, B_ECI, sunPos_ECI, attParams, scParams, epoch)
%==========================================================================
% simulateAttitude: High-fidelity attitude dynamics integration with 
%                   optimized environmental models.
%
% Description:
%   Integrates the spacecraft rotational kinematics (quaternions) and 
%   dynamics (Euler equations) using a 4th-Order Runge-Kutta (RK4) method.
%   Computes external perturbation torques (Gravity Gradient, Aerodynamic 
%   Drag, Solar Radiation Pressure, Magnetic Dipole) and active control 
%   torques (PID Nadir pointing).
%
% Inputs:
%   t          - Time vector                                        [s], 1xN
%   rECI       - Spacecraft position in ECI frame                   [m], 3xN
%   vECI       - Spacecraft velocity in ECI frame                 [m/s], 3xN
%   B_ECI      - Earth magnetic field in ECI frame                 [nT], 3xN
%   sunPos_ECI - Sun position vector in ECI frame                   [m], 3xN
%   epoch      - Simulation start time (datetime object, UTC)
%
%   attParams  - Attitude control parameters structure containing:
%                * .I       : Moment of inertia tensor         [kg·m^2], 3x3
%                * .dipole  : Residual magnetic dipole moment   [A·m^2], 3x1
%                * .Kp      : Proportional control gain             [-], 1x1
%                * .Ki      : Integral control gain (if used)       [-], 1x1
%                * .Kd      : Derivative control gain               [-], 1x1
%
%   scParams   - Spacecraft physical parameters structure containing:
%                * .mass    : Total spacecraft mass                [kg], 1x1
%                * .areaDrag: Cross-sectional area for aero drag  [m^2], 1x1
%                * .areaSRP : Cross-sectional area for SRP        [m^2], 1x1
%                * .Cd      : Aerodynamic drag coefficient          [-], 1x1
%                * .Cr      : Solar radiation pressure coeff.       [-], 1x1
%
% Outputs:
%   qTrue      - True attitude quaternion history (ECI to Body)        , 4xN
%   omegaTrue  - True angular velocity history (body frame)     [rad/s], 3xN
%   torques    - Struct containing torque histories:
%                * .gg      : Gravity gradient torque             [N·m], 3xN
%                * .drag    : Aerodynamic drag torque             [N·m], 3xN
%                * .srp     : Solar radiation pressure torque     [N·m], 3xN
%                * .mag     : Residual magnetic dipole torque     [N·m], 3xN
%                * .ctrl    : Active control torque               [N·m], 3xN
%                * .total   : Total sum of torques                [N·m], 3xN
%==========================================================================

    N = numel(t);
    muEarth = 3.986004418e14; % Earth's gravitational parameter [m^3/s^2]
    REarth  = 6378137.0;      % Earth's equatorial radius       [m]
    wEarth = 7.2921150e-5;    % Earth rotation rate             [rad/s]

    %% ====================================================================
    % 1. PRE-COMPUTE KINEMATICS & ECLIPSE
    % =====================================================================
    % Mean motion (used for LVLH feed-forward term)
    r_norm = sqrt(sum(rECI.^2, 1));
    n      = sqrt(muEarth ./ r_norm.^3);
    
    % LVLH Frame DCM elements (Z-Nadir, Y-CrossTrack, X-AlongTrack)
    yLVLH      = cross(rECI, vECI, 1);
    yLVLH_norm = sqrt(sum(yLVLH.^2, 1));
    yLVLH      = -yLVLH ./ yLVLH_norm;

    zLVLH      = -rECI ./ r_norm;

    xLVLH      = cross(yLVLH, zLVLH, 1);
    
    % Eclipse state (Cylindrical shadow model)
    sunNorm  = sqrt(sum(sunPos_ECI.^2, 1));
    sunDir   = sunPos_ECI ./ sunNorm;
    proj     = sum(rECI .* sunDir, 1);
    perpVec  = rECI - proj .* sunDir;
    perpDist = sqrt(sum(perpVec.^2, 1));
    eclipse  = (proj < 0) & (perpDist < REarth);
    
    % Relative velocity to atmosphere (including Earth rotation wind)
    vAtmECI = [-wEarth * rECI(2,:); 
                wEarth * rECI(1,:); 
                zeros(1, N)];
    vRelECI = vECI - vAtmECI;

    %% ====================================================================
    % 2. PRE-COMPUTE ATMOSPHERIC DENSITY
    % =====================================================================
    % Atmosphere changes slowly. Evaluate NRLMSISE-00 every 10s and interpolate.
    dtEnv    = 10; 
    dsFactor = max(1, round(dtEnv / mean(diff(t))));
    idxEnv   = 1:dsFactor:N;
    
    % Ensure the last point is included to prevent extrapolation errors
    if idxEnv(end) ~= N
        idxEnv = [idxEnv, N]; 
    end
    
    tEnv   = t(idxEnv);
    rhoEnv = zeros(1, length(idxEnv));
    
    % Calculate ECEF DCM only for the downsampled points
    DCM_ECI2ECEF_env = dcmeci2ecef('IAU-2000/2006', epoch + seconds(tEnv'));
    
    for i = 1:length(idxEnv)
        idx = idxEnv(i);
        rECEF = squeeze(DCM_ECI2ECEF_env(:,:,i)) * rECI(:, idx);
        LLA = ecef2lla(rECEF', 'WGS84');
        alt = LLA(3);
        
        if alt > 1000e3
            % Density is negligible above 1000 km altitude
            rhoEnv(i) = 0;
        else
            currTime = epoch + seconds(t(idx));
            doy      = day(currTime, 'dayofyear');
            doySec   = hour(currTime)*3600 + minute(currTime)*60 + second(currTime);
            
            % F10.7=150 (Average solar activity), AP=4 (Quiet geomagnetic activity)
            [~, rhoOut] = atmosnrlmsise00(alt, LLA(1), LLA(2), year(currTime), ...
                                          doy, doySec, 150, 150, 4*ones(1,7), 'Oxygen');
            rhoEnv(i) = rhoOut(6); % Extract total mass density [kg/m^3]
        end
    end
    
    % Interpolate density back to the full high-frequency time array
    rho = interp1(tEnv, rhoEnv, t, 'linear', 'extrap');

    %% ====================================================================
    % 3. ATTITUDE INTEGRATION
    % =====================================================================
    
    % Initial attitude: Align perfectly with LVLH frame at t=0
    DCM_LVLH0 = [xLVLH(:,1)'; yLVLH(:,1)'; zLVLH(:,1)'];
    q0        = dcm2quat(DCM_LVLH0);
    omega0    = [0; 0; 0]; % Start at rest relative to ECI frame
    
    % Pre-allocate outputs
    qTrue     = zeros(4,N);
    omegaTrue = zeros(3,N);
    torques   = struct(  'gg', zeros(3,N),  'drag', zeros(3,N), ...
                        'srp', zeros(3,N),   'mag', zeros(3,N), ...
                       'ctrl', zeros(3,N), 'total', zeros(3,N));
                   
    qTrue(:,1)     = q0;
    omegaTrue(:,1) = omega0;
    
    % High-frequency integration loop
    for k = 1:(N-1)
        dt = t(k+1) - t(k);
        
        % Fast lookup for current state environment
        env.rECI     = rECI(:,k);
        env.B_ECI    = B_ECI(:,k);
        env.sunECI   = sunPos_ECI(:,k);
        env.vRelECI  = vRelECI(:,k);
        env.rho      = rho(k);
        env.eclipse  = eclipse(k);
        env.n        = n(k);
        env.DCM_LVLH = [xLVLH(:,k)'; yLVLH(:,k)'; zLVLH(:,k)'];
        
        q_k = qTrue(:,k);
        w_k = omegaTrue(:,k);
        
        % RK4 Core (assuming environment is constant over dt for hyper-speed)
        [k1_q, k1_w, tau_k] = computeAttitudeDerivatives(q_k, w_k, env, attParams, scParams);
        
        % Store torques for the current step
        torques.gg(:,k)    = tau_k.gg;
        torques.drag(:,k)  = tau_k.drag;
        torques.srp(:,k)   = tau_k.srp;
        torques.mag(:,k)   = tau_k.mag;
        torques.ctrl(:,k)  = tau_k.ctrl;
        torques.total(:,k) = tau_k.total;
        
        % RK4 intermediate steps
        qMid = q_k + 1/2*dt*k1_q;
        qMid = qMid / norm(qMid);
        wMid = w_k + 1/2*dt*k1_w;
        [k2_q, k2_w, ~] = computeAttitudeDerivatives(qMid, wMid, env, attParams, scParams);
        
        qMid2 = q_k + 1/2*dt*k2_q;
        qMid2 = qMid2 / norm(qMid2);
        wMid2 = w_k + 1/2*dt*k2_w;
        [k3_q, k3_w, ~] = computeAttitudeDerivatives(qMid2, wMid2, env, attParams, scParams);
        
        qEnd = q_k + dt*k3_q;
        qEnd = qEnd / norm(qEnd);
        wEnd = w_k + dt*k3_w;
        [k4_q, k4_w, ~] = computeAttitudeDerivatives(qEnd, wEnd, env, attParams, scParams);
        
        % Update state
        qTrue(:,k+1)     = q_k + 1/6*dt*(k1_q + 2*k2_q + 2*k3_q + k4_q);
        qTrue(:,k+1)     = qTrue(:,k+1) / norm(qTrue(:,k+1));
        omegaTrue(:,k+1) = w_k + 1/6*dt*(k1_w + 2*k2_w + 2*k3_w + k4_w);
    end
    
    % Evaluate torque at the final step to complete the arrays
    env.rECI     = rECI(:,N);
    env.B_ECI    = B_ECI(:,N);
    env.sunECI   = sunPos_ECI(:,N);
    env.vRelECI  = vRelECI(:,N);
    env.rho      = rho(N);
    env.eclipse  = eclipse(N);
    env.n        = n(N);
    env.DCM_LVLH = [xLVLH(:,N)'; yLVLH(:,N)'; zLVLH(:,N)'];
    
    [~, ~, tauFinal] = computeAttitudeDerivatives(qTrue(:,N), omegaTrue(:,N), env, attParams, scParams);
    torques.gg(:,N)    = tauFinal.gg;
    torques.drag(:,N)  = tauFinal.drag;
    torques.srp(:,N)   = tauFinal.srp;
    torques.mag(:,N)   = tauFinal.mag;
    torques.ctrl(:,N)  = tauFinal.ctrl;
    torques.total(:,N) = tauFinal.total;
end

%% =======================================================================
% LOCAL HELPER: Attitude Derivatives & Torques Computation
% ========================================================================
function [qDot, omegaDot, torques] = computeAttitudeDerivatives(q, omega, env, attParams, scParams)
%==========================================================================
% computeAttitudeDerivatives: Calculates rotational kinematics and dynamics
%                             for a single RK4 integration step.
%
% Description:
%   Evaluates the time derivatives of the attitude quaternion (kinematics) 
%   and angular velocity vector (Euler's equations). Computes instantaneous 
%   environmental disturbance torques (Gravity Gradient, Aerodynamic Drag, 
%   Solar Radiation Pressure, Magnetic Dipole) and the active PD control 
%   torque required for Nadir tracking relative to the LVLH frame.
%
% Inputs:
%   q          - Current attitude quaternion (ECI to Body)             , 4x1
%   omega      - Current angular velocity in body frame         [rad/s], 3x1
%   attParams  - Struct with attitude control parameters (see above)
%   scParams   - Struct with spacecraft physical parameters (see above)
%
%   env        - Struct with instantaneous environmental conditions:
%                * .rECI    : Position in ECI frame                 [m], 3x1
%                * .B_ECI   : Earth magnetic field in ECI          [nT], 3x1
%                * .sunECI  : Sun position vector in ECI            [m], 3x1
%                * .vRelECI : Velocity relative to atmosphere     [m/s], 3x1
%                * .rho     : Atmospheric density              [kg/m^3], 1x1
%                * .eclipse : Eclipse flag (1: shadow, 0: sun)      [-], 1x1
%                * .n       : Mean motion                       [rad/s], 1x1
%                * .DCM_LVLH: LVLH frame Direction Cosine Matrix    [-], 3x3
%
% Outputs:
%   qDot       - Time derivative of the attitude quaternion       [1/s], 4x1
%   omegaDot   - Angular acceleration in the body frame       [rad/s^2], 3x1
%   torques    - Struct containing instantaneous torques          [N·m]
%==========================================================================

    muEarth = 3.986004418e14; % [m^3/s^2]
    
    % Current Direction Cosine Matrix (ECI to Body)
    DCM_ECI2B = quat2dcm(q);
    
    % --- 1. GRAVITY GRADIENT TORQUE ---
    rBody = DCM_ECI2B * env.rECI;
    rNorm = norm(rBody);
    rDir  = rBody / rNorm;
    tauGG = (3*muEarth / rNorm^3) * cross(rDir, attParams.I * rDir);
    
    % --- 2. AERODYNAMIC DRAG TORQUE ---
    if env.rho > 0
        vRel_body = DCM_ECI2B * env.vRelECI;
        vNorm    = norm(vRel_body);
        if vNorm > 1e-3
            % Assume Center of Pressure is offset +5 cm along X-axis
            rCP = [0.05; 0; 0]; 
            tauDrag = 0.5 * env.rho * vNorm^2 * scParams.Cd * scParams.areaDrag * cross(rCP, -(vRel_body/vNorm));
        else
            tauDrag = zeros(3,1);
        end
    else
        tauDrag = zeros(3,1);
    end
    
    % --- 3. SOLAR RADIATION PRESSURE TORQUE ---
    if ~env.eclipse
        sunBody = DCM_ECI2B * env.sunECI;
        sunNorm = norm(sunBody);
        if sunNorm > 1e-3
            P_sun = 4.56e-6 * (1.496e11 / sunNorm)^2; % Solar pressure [N/m²]
            % Assume Center of Pressure is offset +10 cm along Y-axis
            rCP = [0; 0.1; 0]; 
            tauSRP = P_sun * scParams.Cr * scParams.areaSRP * cross(rCP, sunBody/sunNorm);
        else
            tauSRP = zeros(3,1);
        end
    else
        tauSRP = zeros(3,1);
    end
    
    % --- 4. MAGNETIC DIPOLE TORQUE ---
    % Convert magnetic field from nT to Tesla
    B_body = DCM_ECI2B * env.B_ECI * 1e-9;
    tauMAG = cross(attParams.dipole, B_body);
    
    % --- 5. CONTROL TORQUE (Nadir Tracking) ---
    % Determine orientation error relative to LVLH frame
    DCMErr   = DCM_ECI2B * env.DCM_LVLH';
    thetaErr = [DCMErr(3,2) - DCMErr(2,3); 
                DCMErr(1,3) - DCMErr(3,1); 
                DCMErr(2,1) - DCMErr(1,2)] / 2;
                 
    % Determine angular velocity error relative to LVLH frame
    omega_LVLH_body = DCM_ECI2B * env.DCM_LVLH' * [0; -env.n; 0];
    omegaErr        = omega - omega_LVLH_body;
    
    % PD Control Law
    tauCtrl = -attParams.Kp * thetaErr - attParams.Kd * omegaErr;
    
    % --- 6. TOTAL TORQUE & KINEMATICS ---
    tauTotal = tauGG + tauDrag + tauSRP + tauMAG + tauCtrl;
    
    torques.gg    = tauGG;
    torques.drag  = tauDrag;
    torques.srp   = tauSRP;
    torques.mag   = tauMAG;
    torques.ctrl  = tauCtrl;
    torques.total = tauTotal;
    
    % Attitude Kinematics (Quaternion Derivative)
    qDot = 1/2 * skewMatrix(omega) * q;
    
    % Rigid Body Dynamics (Euler's Equation)
    omegaDot = attParams.I \ (tauTotal - cross(omega, attParams.I * omega));
end