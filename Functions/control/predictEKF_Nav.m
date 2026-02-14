function EKF = predictEKF_Nav(EKF, forceMeas, q_ECI2B, dt)
%==========================================================================
% predictEKF_Nav: Propagates the PV-EKF state and covariance through the 
%                 prediction step using accelerometer measurements and 
%                 orbital dynamics.
%
% Inputs:
%    EKF       - PV-EKF structure with fields:
%        .x    - State vector (9x1): [r_ECI; v_ECI; biasAccel]
%        .P    - State covariance matrix (9x9)
%        .Q    - Discrete process noise covariance matrix (9x9)
%
%    forceMeas - Measured specific force from accelerometer [m/s^2] (3x1)
%                (This is in the Body frame)
%    q_ECI2B   - Current attitude quaternion (ECI to Body) (4x1)
%    dt        - Time step [s]
%
% Outputs:
%    EKF       - Updated EKF structure (propagated state and covariance).
%
% Algorithm:
%    1. State Propagation: Integrates position and velocity using central
%       gravity + J2 perturbation, plus the bias-corrected specific force
%       rotated to the ECI frame.
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
    ba = EKF.x(7:9);
    
    r_norm = norm(r);

    %% ====================================================================
    % 1. ACCELERATION COMPUTATION
    % =====================================================================
    
    % --- A. Gravitational Acceleration (Central Body + J2) ---
    a_g = -muEarth / r_norm^3 * r;
    
    % J2 Perturbation
    z2 = (r(3) / r_norm)^2;
    J2_fac = 1.5 * J2Earth * muEarth * REarth^2 / r_norm^5;
    a_J2 = zeros(3,1);
    a_J2(1) = J2_fac * r(1) * (5*z2 - 1);
    a_J2(2) = J2_fac * r(2) * (5*z2 - 1);
    a_J2(3) = J2_fac * r(3) * (5*z2 - 3);
    
    g_total = a_g + a_J2;
    
    % --- B. Specific Force (Accelerometer) ---
    % Correct measurement with current estimated bias
    f_body = forceMeas - ba;
    
    % Rotate specific force from Body to ECI frame
    DCM_ECI2B = quat2dcm(q_ECI2B);
    DCM_B2ECI = DCM_ECI2B';
    f_ECI     = DCM_B2ECI * f_body;
    
    % --- C. Total Acceleration ---
    a_total = g_total + f_ECI;

    %% ====================================================================
    % 2. STATE PROPAGATION (Numerical Integration)
    % =====================================================================
    
    % 2nd-order Euler integration for position
    EKF.x(1:3) = r + v * dt + 0.5 * a_total * dt^2;
    
    % 1st-order Euler integration for velocity
    EKF.x(4:6) = v + a_total * dt;
    
    % Bias is modeled as a random walk, so its deterministic derivative is 0
    % EKF.x(7:9) remains unchanged in the prediction step

    %% ====================================================================
    % 3. COVARIANCE PROPAGATION
    % =====================================================================
    
    I3 = eye(3);
    O3 = zeros(3,3);
    
    % Gravity Gradient Tensor (Partial derivative of gravity w.r.t position)
    % (J2 gradient is negligible for covariance propagation, central body is enough)
    G = (muEarth / r_norm^5) * (3 * (r * r') - (r_norm^2) * I3);
    
    % Continuous-time Jacobian Matrix (F)
    % State: [r; v; ba]
    % d(r_dot)/dr = 0, d(r_dot)/dv = I, d(r_dot)/dba = 0
    % d(v_dot)/dr = G, d(v_dot)/dv = 0, d(v_dot)/dba = -DCM_B2ECI
    % d(ba_dot)/dr = 0, d(ba_dot)/dv = 0, d(ba_dot)/dba = 0
    
    F = [ O3,  I3,        O3;
           G,  O3, -DCM_B2ECI;
          O3,  O3,        O3 ];
      
    % Discrete-time State Transition Matrix (STM) - 1st order Taylor expansion
    Phi = eye(9) + F * dt;
    
    % Propagate Covariance (P = Phi * P * Phi' + Q)
    % Note: EKF.Q was already discretized in initializeEKF_Nav.m, so no '*dt' here
    EKF.P = Phi * EKF.P * Phi' + EKF.Q;
    
    % Force symmetry to prevent numerical drift over time
    EKF.P = (EKF.P + EKF.P') / 2;

end