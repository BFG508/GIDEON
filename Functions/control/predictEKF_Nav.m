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
    bA = EKF.x(7:9);
    
    rNorm = norm(r);

    %% ====================================================================
    % 1. ACCELERATION COMPUTATION
    % =====================================================================
    
    % --- A. Gravitational Acceleration (Central Body + J2) ---
    % Acceleration due to central body
    aG = -muEarth * r / rNorm^3;
    
    % Acceleration due to J2 (Oblateness)
    % Factor common to all components
    J2Const  = (1.5 * J2Earth * muEarth * REarth^2) / rNorm^5;
    zSquared = (r(3) / rNorm)^2;
    
    aJ2   = zeros(3,1);
    aJ2(1) = J2Const * r(1) * (5 * zSquared - 1); % x-component
    aJ2(2) = J2Const * r(2) * (5 * zSquared - 1); % y-component
    aJ2(3) = J2Const * r(3) * (5 * zSquared - 3); % z-component
    
    gTotal = aG + aJ2;
    
    % --- B. Specific Force (Accelerometer) ---
    % Correct measurement with current estimated bias
    f_body = forceMeas - bA;
    
    % Rotate specific force from Body to ECI frame
    DCM_ECI2B = quat2dcm(q_ECI2B);
    DCM_B2ECI = DCM_ECI2B';
    f_ECI     = DCM_B2ECI * f_body;
    
    % --- C. Total Acceleration ---
    aTotal = gTotal + f_ECI;

    %% ====================================================================
    % 2. STATE PROPAGATION (Numerical Integration)
    % =====================================================================
    
    % 2nd-order Euler integration for position
    EKF.x(1:3) = r + v * dt + 1/2 * aTotal * dt^2;
    
    % 1st-order Euler integration for velocity
    EKF.x(4:6) = v + aTotal * dt;
    
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