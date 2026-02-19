function EKF = predictEKF_Nav(EKF, forceMeas, q_ECI2B, dt)
%==========================================================================
% predictEKF_Nav: Propagates the PV-EKF state and covariance through the 
%                 prediction step using accelerometer measurements and 
%                 orbital dynamics.
%
% Inputs:
%    EKF       - PV-EKF structure with fields:
%        .x    - State vector (9x1): [rECI; vECI; biasAccel]
%        .P    - State covariance matrix (9x9)
%        .Q    - Discrete process noise covariance matrix (9x9)
%    forceMeas - Measured specific force from accelerometer [m/s²] (3x1)
%    q_ECI2B   - Current attitude quaternion (ECI to Body) (4x1)
%    dt        - Time step [s]
%
% Outputs:
%    EKF       - Updated EKF structure propagated.
%
% Algorithm:
%    1. State Propagation: Integrates position and velocity using RK4 with
%       central gravity + J2 perturbation, plus the bias-corrected specific
%       force rotated to the ECI frame.
%    2. Covariance Propagation: Evaluates the continuous-time Jacobian (F)
%       including the gravity gradient tensor, computes the State Transition
%       Matrix (Phi), and propagates P.
%==========================================================================

    % Constants
    muEarth = 3.986004418e14;  % [m^3/s^2]
    REarth  = 6378137.0;       % [m]
    J2Earth = 1.08262668e-3;   % [-]

    % Extract current state
    r  = EKF.x(1:3);
    v  = EKF.x(4:6);
    bA = EKF.x(7:9);
    
    %% ====================================================================
    % 1. SPECIFIC FORCE COMPUTATION
    % =====================================================================
    
    % Correct measurement with current estimated bias
    f_body = forceMeas - bA;
    
    % Rotate specific force from Body to ECI frame
    DCM_ECI2B = quat2dcm(q_ECI2B);
    DCM_B2ECI = DCM_ECI2B';
    f_ECI     = DCM_B2ECI * f_body;
    
    %% ====================================================================
    % 2. STATE PROPAGATION
    % =====================================================================
    
    % Create a fast anonymous function to evaluate gravity 
    % (Central Body + J2) at any given position "p"
    calcGrav = @(p) (-muEarth / norm(p)^3) * p + ...
        ((3/2 * J2Earth * muEarth * REarth^2) / norm(p)^5) * ...
        [p(1) * (5 * (p(3) / norm(p))^2 - 1);
         p(2) * (5 * (p(3) / norm(p))^2 - 1);
         p(3) * (5 * (p(3) / norm(p))^2 - 3)]; 
         
    % Note: We assume the specific force (f_ECI) measured by the IMU 
    % remains constant during this small time step (dt)
    
    % --- Step k1 ---
    v1 = v;
    a1 = calcGrav(r) + f_ECI;
    
    % --- Step k2 ---
    r2 = r + 1/2 * v1 * dt;
    v2 = v + 1/2 * a1 * dt;
    a2 = calcGrav(r2) + f_ECI;
    
    % --- Step k3 ---
    r3 = r + 1/2 * v2 * dt;
    v3 = v + 1/2 * a2 * dt;
    a3 = calcGrav(r3) + f_ECI;
    
    % --- Step k4 ---
    r4 = r + v3 * dt;
    v4 = v + a3 * dt;
    a4 = calcGrav(r4) + f_ECI;
    
    % --- Final EKF state update ---
    EKF.x(1:3) = r + (dt / 6) * (v1 + 2*v2 + 2*v3 + v4);
    EKF.x(4:6) = v + (dt / 6) * (a1 + 2*a2 + 2*a3 + a4);
    
    % Bias is modeled as a random walk, so its deterministic derivative is 0
    % EKF.x(7:9) remains unchanged in the prediction step
    
    %% ====================================================================
    % 3. COVARIANCE PROPAGATION
    % =====================================================================
   
    % Discrete-time State Transition Matrix (STM) - 1st order Taylor expansion
    Phi = computeSTM_Nav(r, q_ECI2B, dt);
    
    % Propagate Covariance (P = Phi * P * Phi' + Q)
    EKF.P = Phi * EKF.P * Phi' + EKF.Q;
    
    % Force symmetry to prevent numerical drift over time
    EKF.P = (EKF.P + EKF.P') / 2;

end