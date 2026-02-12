function truth = generateTruthTrajectory(t, epoch, orbitalElements, scParams, attParams)
%==========================================================================
% generateTruthTrajectory - Generate realistic spacecraft trajectory with
%                           attitude dynamics and environmental perturbations
%
% INPUTS:
%   t               - Time vector [s], Nx1 or 1xN
%   epoch           - Simulation start time (datetime object, UTC)
%   orbitalElements - Orbital parameters (see simulateOrbit.m)
%   scParams        - Spacecraft parameters (mass, areas, coefficients)
%   attParams       - Attitude parameters (struct):
%       .I          - Inertia tensor [kg·m²], 3x3
%       .dipole     - Magnetic dipole moment [A·m²], 3x1 (body frame)
%       .Kp         - Proportional gain for attitude control [N·m]
%       .Kd         - Derivative gain for attitude control [N·m·s]
%       .mode       - 'LVLH' (Nadir Pointing) or 'INERTIAL'
%
% OUTPUTS:
%   truth - Structure with fields:
%       .t          - Time vector [s], 1xN
%       .rECI       - Position in ECI [m], 3xN
%       .vECI       - Velocity in ECI [m/s], 3xN
%       .qTrue      - Attitude quaternion (ECI to Body), 4xN
%       .omegaTrue  - Angular velocity (body frame) [rad/s], 3xN
%       .B_ECI      - Magnetic field in ECI [nT], 3xN
%       .torques    - Applied torques (body frame) [N·m], struct
%==========================================================================

    fprintf('\n=== Generating Truth Trajectory ===\n');
    
    t = t(:)';  % Ensure row vector
    N = numel(t);
    
    %% ===================================================================
    % 1. ORBITAL PROPAGATION (Position and Velocity)
    % ====================================================================
    fprintf(' Propagating orbit ... ');
    [rECI, vECI] = simulateOrbit(t, orbitalElements, scParams);
    fprintf('Done.\n');
    
    %% ===================================================================
    % 2. EARTH MAGNETIC FIELD (IGRF-13) with ECI↔ECEF Downsampling
    % ====================================================================
    fprintf(' Computing geomagnetic field (IGRF-13) ... \n');
    
    % Downsample DCM Calculation (ECI → ECEF transformation is slow)
    dt_original = mean(diff(t));
    dt_dcm      = 1.0;  % Calculate DCM every 1 second (instead of IMU rate)
    downsample_factor = round(dt_dcm / dt_original);
    
    % Downsampled time vector
    t_dcm = t(1:downsample_factor:end);
    N_dcm = length(t_dcm);
    
    fprintf('   Original samples: %d at %.2f Hz\n', N, 1/dt_original);
    fprintf('   DCM samples:      %d at %.2f Hz (downsample factor: %d)\n', ...
            N_dcm, 1/dt_dcm, downsample_factor);
    
    % Compute DCM_ECI2ECEF at downsampled rate
    fprintf('   Computing ECI→ECEF transformation ... ');
    tic;
    DCM_ECI2ECEF_downsampled = dcmeci2ecef('IAU-2000/2006', epoch + seconds(t_dcm'));
    t_dcm_calc = toc;
    fprintf('Done (%.2f s)\n', t_dcm_calc);
    
    % Interpolate DCM to original time vector
    fprintf('   Interpolating DCM to full time vector ... ');
    tic;
    DCM_ECI2ECEF = zeros(3, 3, N);
    
    for i = 1:3
        for j = 1:3
            % Extract time series for this DCM element
            dcm_ij = squeeze(DCM_ECI2ECEF_downsampled(i,j,:));
            
            % Interpolate to full time vector (linear is sufficient)
            DCM_ECI2ECEF(i,j,:) = interp1(t_dcm, dcm_ij, t, 'linear', 'extrap');
        end
    end
    t_interp = toc;
    fprintf('Done (%.2f s)\n', t_interp);
    
    % Convert ECI positions to ECEF using interpolated DCM
    fprintf('   Transforming positions ECI→ECEF ... ');
    rECEF = zeros(3, N);
    for k = 1:N
        rECEF(:,k) = squeeze(DCM_ECI2ECEF(:,:,k)) * rECI(:,k);
    end
    fprintf('Done.\n');
    
    % Convert ECEF to LLA (Latitude, Longitude, Altitude)
    fprintf('   Converting ECEF→LLA ... ');
    LLA = ecef2lla(rECEF', 'WGS84');
    lat = LLA(:,1)';  % [deg]
    lon = LLA(:,2)';  % [deg]
    alt = LLA(:,3)';  % [m]
    fprintf('Done.\n');
    
    % Compute Magnetic Field using IGRF-13
    fprintf('   Evaluating IGRF-13 model ... ');
    XYZ_NED   = igrfmagm(alt/1e3, lat, lon, 2025.0 * ones(1,N), 13);  % alt in km
    fprintf('Done.\n');
    
    % Convert NED → ECEF → ECI
    fprintf('   Transforming magnetic field NED→ECEF→ECI ... ');
    [B_ECEFx, B_ECEFy, B_ECEFz] = ned2ecefv(XYZ_NED(:,1), XYZ_NED(:,2), XYZ_NED(:,3), ...
                                            lat(:), lon(:));
    B_ECEF = [B_ECEFx'; B_ECEFy'; B_ECEFz'];
    
    % Rotate ECEF → ECI using transposed DCM
    B_ECI = zeros(3, N);
    for k = 1:N
        DCM_ECEF2ECI = squeeze(DCM_ECI2ECEF(:,:,k))';  % Transpose
        B_ECI(:,k) = DCM_ECEF2ECI * B_ECEF(:,k);
    end
    fprintf('Done.\n');
    
    fprintf(' Magnetic field computation complete.\n');
    
    %% ===================================================================
    % 3. SUN POSITION (for SRP torque)
    % ====================================================================
    fprintf(' Computing Sun ephemeris ... ');
    sunPos_ECI = zeros(3,N);
    for k = 1:N
        sunPos_ECI(:,k) = computeSunPosition(epoch + seconds(t(k)));
    end
    fprintf('Done.\n');
    
    %% ===================================================================
    % 4. ATTITUDE DYNAMICS INTEGRATION
    % ====================================================================
    fprintf(' Integrating attitude dynamics (Nadir Pointing + perturbations) ... ');
    
    % Initial attitude: Align with LVLH frame
    DCM_LVLH0 = computeLVLH_DCM(rECI(:,1), vECI(:,1));
    q0        = dcm2quat(DCM_LVLH0);
    omega0    = [0; 0; 0];  % Start at rest (relative to LVLH)
    
    % Pre-allocate outputs
    qTrue     = zeros(4,N);
    omegaTrue = zeros(3,N);
    torques   = struct('gg', zeros(3,N), 'drag', zeros(3,N), ...
                       'srp', zeros(3,N), 'mag', zeros(3,N), ...
                       'control', zeros(3,N), 'total', zeros(3,N));
    
    % Initial state
    qTrue(:,1)     = q0;
    omegaTrue(:,1) = omega0;
    
    % Integrate using RK4 (explicit, simple)
    for k = 1:(N-1)
        dt = t(k+1) - t(k);
        
        % Current state
        q_k     = qTrue(:,k);
        omega_k = omegaTrue(:,k);
        r_k     = rECI(:,k);
        v_k     = vECI(:,k);
        B_k     = B_ECI(:,k);
        sun_k   = sunPos_ECI(:,k);
        
        % Compute torques and derivatives
        [qDot_k, omegaDot_k, tau_k] = attitudeDynamics(q_k, omega_k, r_k, v_k, B_k, sun_k, ...
                                                        attParams, scParams);
        
        % Store torques
        torques.gg(:,k)      = tau_k.gg;
        torques.drag(:,k)    = tau_k.drag;
        torques.srp(:,k)     = tau_k.srp;
        torques.mag(:,k)     = tau_k.mag;
        torques.control(:,k) = tau_k.control;
        torques.total(:,k)   = tau_k.total;
        
        % RK4 integration
        k1_q = qDot_k;
        k1_w = omegaDot_k;
        
        q_mid = q_k + 0.5*dt*k1_q;
        q_mid = q_mid / norm(q_mid);
        omega_mid = omega_k + 0.5*dt*k1_w;
        [qDot_mid, omegaDot_mid, ~] = attitudeDynamics(q_mid, omega_mid, r_k, v_k, B_k, sun_k, ...
                                                        attParams, scParams);
        k2_q = qDot_mid;
        k2_w = omegaDot_mid;
        
        q_mid2 = q_k + 0.5*dt*k2_q;
        q_mid2 = q_mid2 / norm(q_mid2);
        omega_mid2 = omega_k + 0.5*dt*k2_w;
        [qDot_mid2, omegaDot_mid2, ~] = attitudeDynamics(q_mid2, omega_mid2, r_k, v_k, B_k, sun_k, ...
                                                          attParams, scParams);
        k3_q = qDot_mid2;
        k3_w = omegaDot_mid2;
        
        q_end = q_k + dt*k3_q;
        q_end = q_end / norm(q_end);
        omega_end = omega_k + dt*k3_w;
        [qDot_end, omegaDot_end, ~] = attitudeDynamics(q_end, omega_end, r_k, v_k, B_k, sun_k, ...
                                                        attParams, scParams);
        k4_q = qDot_end;
        k4_w = omegaDot_end;
        
        % Update state
        qTrue(:,k+1)     = q_k + (dt/6)*(k1_q + 2*k2_q + 2*k3_q + k4_q);
        qTrue(:,k+1)     = qTrue(:,k+1) / norm(qTrue(:,k+1));
        omegaTrue(:,k+1) = omega_k + (dt/6)*(k1_w + 2*k2_w + 2*k3_w + k4_w);
    end
    
    % Final torque evaluation
    [~, ~, tau_final] = attitudeDynamics(qTrue(:,N), omegaTrue(:,N), rECI(:,N), vECI(:,N), ...
                                          B_ECI(:,N), sunPos_ECI(:,N), attParams, scParams);
    torques.gg(:,N)      = tau_final.gg;
    torques.drag(:,N)    = tau_final.drag;
    torques.srp(:,N)     = tau_final.srp;
    torques.mag(:,N)     = tau_final.mag;
    torques.control(:,N) = tau_final.control;
    torques.total(:,N)   = tau_final.total;
    
    fprintf('Done.\n');
    
    %% ===================================================================
    % 5. ASSEMBLE OUTPUT STRUCTURE
    % ====================================================================
    truth.t         = t;
    truth.rECI      = rECI;
    truth.vECI      = vECI;
    truth.qTrue     = qTrue;
    truth.omegaTrue = omegaTrue;
    truth.B_ECI     = B_ECI;
    truth.torques   = torques;
    
    fprintf('=== Truth Trajectory Complete ===\n');
    fprintf(' Mean altitude:        %.1f km\n', mean(vecnorm(rECI))/1e3 - 6378);
    fprintf(' Mean angular rate:    [%.4f, %.4f, %.4f] deg/s\n', ...
            mean(rad2deg(abs(omegaTrue)),2));
    fprintf(' Mean magnetic field:  %.1f nT\n', mean(vecnorm(B_ECI)));
end

%% =======================================================================
% HELPER FUNCTION: Attitude Dynamics with Perturbation Torques
% ========================================================================
function [qDot, omegaDot, torques] = attitudeDynamics(q, omega, r, v, B, sunPos, attParams, scParams)
    % Constants
    muEarth = 3.986004418e14;  % [m^3/s^2]
    
    % Extract parameters
    I      = attParams.I;
    dipole = attParams.dipole;
    Kp     = attParams.Kp;
    Kd     = attParams.Kd;
    
    % Current DCM (ECI to Body)
    DCM_ECI2B = quat2dcm(q);
    
    %% 1. GRAVITY GRADIENT TORQUE
    r_body = DCM_ECI2B * r;
    r_norm = norm(r);
    r_hat  = r_body / norm(r_body);
    
    tau_gg = (3*muEarth / r_norm^3) * cross(r_hat, I * r_hat);
    
    %% 2. AERODYNAMIC DRAG TORQUE (simplified)
    % Assume center of pressure offset from center of mass
    rho = computeAtmosphericDensity(r_norm - 6378137);  % Exponential model
    v_body = DCM_ECI2B * v;
    v_norm = norm(v_body);
    
    if v_norm > 1e-3
        v_hat = v_body / v_norm;
        r_cp_cm = [0.05; 0; 0];  % 5 cm offset in +X direction
        tau_drag = 0.5 * rho * v_norm^2 * scParams.Cd * scParams.areaDrag * cross(r_cp_cm, -v_hat);
    else
        tau_drag = zeros(3,1);
    end
    
    %% 3. SOLAR RADIATION PRESSURE TORQUE
    sun_body = DCM_ECI2B * sunPos;
    sun_norm = norm(sun_body);
    
    % Check if in eclipse
    eclipse = checkEclipse(r, sunPos);
    
    if sun_norm > 1e-3 && ~eclipse
        sun_hat = sun_body / sun_norm;
        P_sun   = 4.56e-6 * (1.496e11 / sun_norm)^2;  % [N/m²]
        r_cp_cm = [0; 0.1; 0];  % 10 cm offset in +Y direction
        tau_srp = P_sun * scParams.Cr * scParams.areaSRP * cross(r_cp_cm, sun_hat);
    else
        tau_srp = zeros(3,1);
    end
    
    %% 4. MAGNETIC DIPOLE TORQUE
    B_body = DCM_ECI2B * B * 1e-9;  % Convert nT to T
    tau_mag = cross(dipole, B_body);
    
    %% 5. CONTROL TORQUE (PD Controller for LVLH Tracking)
    % Desired attitude: LVLH frame
    DCM_LVLH = computeLVLH_DCM(r, v);
    DCM_err  = DCM_ECI2B * DCM_LVLH';  % Error DCM (LVLH to Body)
    
    % Small-angle approximation for error vector
    theta_err = [DCM_err(3,2) - DCM_err(2,3);
                 DCM_err(1,3) - DCM_err(3,1);
                 DCM_err(2,1) - DCM_err(1,2)] / 2;
    
    % Desired angular velocity (orbital rate in LVLH frame)
    r_norm = norm(r);
    n      = sqrt(muEarth / r_norm^3);  % Mean motion
    omega_LVLH_inertial = [0; -n; 0];
    omega_LVLH_body     = DCM_ECI2B * DCM_LVLH' * omega_LVLH_inertial;
    omega_err           = omega - omega_LVLH_body;
    
    % PD control law
    tau_control = -Kp * theta_err - Kd * omega_err;
    
    %% 6. TOTAL TORQUE
    tau_total = tau_gg + tau_drag + tau_srp + tau_mag + tau_control;
    
    % Store individual torques
    torques.gg      = tau_gg;
    torques.drag    = tau_drag;
    torques.srp     = tau_srp;
    torques.mag     = tau_mag;
    torques.control = tau_control;
    torques.total   = tau_total;
    
    %% 7. EQUATIONS OF MOTION
    % Quaternion kinematics
    qDot = 0.5 * skewMatrix(omega) * q;
    
    % Euler's rotational equation: I*ωdot + ω × (I*ω) = τ
    omegaDot = I \ (tau_total - cross(omega, I * omega));
end

%% =======================================================================
% HELPER FUNCTION: Compute LVLH Frame DCM
% ========================================================================
function DCM_ECI2LVLH = computeLVLH_DCM(r, v)
    % LVLH Frame Definition:
    %   x_LVLH: Along velocity direction (ram)
    %   z_LVLH: Nadir (towards Earth center, -r direction)
    %   y_LVLH: Cross-track (completes right-handed system)
    
    z_LVLH = -r / norm(r);           % Nadir
    y_LVLH = -cross(r, v);           % Negative orbit normal
    y_LVLH = y_LVLH / norm(y_LVLH);
    x_LVLH = cross(y_LVLH, z_LVLH);  % Ram direction
    
    DCM_ECI2LVLH = [x_LVLH'; y_LVLH'; z_LVLH'];
end

%% =======================================================================
% HELPER FUNCTION: Sun Position (simplified)
% ========================================================================
function sunPos = computeSunPosition(epoch_datetime)
    % Simplified Sun position (low-fidelity, adequate for torque estimation)
    JD = juliandate(epoch_datetime);
    T  = (JD - 2451545.0) / 36525;  % Julian centuries since J2000
    
    % Mean longitude
    L = 280.460 + 36000.771*T;
    L = mod(L, 360);
    
    % Mean anomaly
    g = 357.528 + 35999.050*T;
    g = mod(g, 360);
    g_rad = deg2rad(g);
    
    % Ecliptic longitude
    lambda = L + 1.915*sin(g_rad) + 0.020*sin(2*g_rad);
    lambda_rad = deg2rad(lambda);
    
    % Obliquity of ecliptic
    epsilon = 23.439 - 0.013*T;
    epsilon_rad = deg2rad(epsilon);
    
    % Distance (AU)
    r_AU = 1.00014 - 0.01671*cos(g_rad) - 0.00014*cos(2*g_rad);
    r_m  = r_AU * 1.496e11;  % Convert to meters
    
    % Sun position in ECI
    sunPos = r_m * [cos(lambda_rad);
                    cos(epsilon_rad)*sin(lambda_rad);
                    sin(epsilon_rad)*sin(lambda_rad)];
end

%% =======================================================================
% HELPER FUNCTION: Eclipse Check
% ========================================================================
function inEclipse = checkEclipse(r_sat, r_sun)
    % Simple cylindrical shadow model
    REarth = 6378137;
    
    % Project satellite position onto Sun-Earth line
    r_sun_hat = r_sun / norm(r_sun);
    proj = dot(r_sat, r_sun_hat);
    
    % If behind Earth and within shadow cylinder
    if proj < 0
        perp_dist = norm(r_sat - proj * r_sun_hat);
        inEclipse = (perp_dist < REarth);
    else
        inEclipse = false;
    end
end

%% =======================================================================
% HELPER FUNCTION: Atmospheric Density (Exponential Model)
% ========================================================================
function rho = computeAtmosphericDensity(altitude)
    % altitude in meters above sea level
    % Simple exponential model
    h_km = altitude / 1000;
    
    if h_km < 200
        rho0 = 2.789e-10;
        H    = 58.515;
    elseif h_km < 300
        rho0 = 1.774e-11;
        H    = 60.828;
    elseif h_km < 500
        rho0 = 3.899e-12;
        H    = 63.822;
    else
        rho0 = 3.561e-13;
        H    = 71.835;
    end
    
    rho = rho0 * exp(-(h_km - 200) / H);
end
