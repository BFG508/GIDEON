function Phi = computeSTM_Nav(rECI, q_ECI2B, dt)
%==========================================================================
% computeSTM_Nav: Computes State Transition Matrix (STM) for PV-EKF error state
%                 propagation using first-order discrete approximation.
%
% Inputs:
%    rECI    - Position vector in ECI [m] (3x1)
%    q_ECI2B - Attitude Quaternion (ECI -> Body) (4x1)
%    dt      - Time step [s]
%
% Outputs:
%    Phi     - Discrete State Transition Matrix (9x9)
%==========================================================================
    
    % Constants
    muEarth = 3.986004418e14; 
    
    r_norm = norm(rECI);
    I3 = eye(3);
    O3 = zeros(3,3);

    % --- 1. Gravity Gradient Tensor (G) ---
    % Partial derivative of gravity accel w.r.t position
    % (J2 gradient is negligible for covariance propagation, central body is enough)
    G = (muEarth / r_norm^5) * (3 * (rECI * rECI') - (r_norm^2) * I3);
    
    % --- 2. Rotation Matrix (Body to ECI) ---
    % Needed to project accelerometer bias uncertainty into ECI frame
    DCM_ECI2B = quat2dcm(q_ECI2B);
    DCM_B2ECI = DCM_ECI2B'; 

    % --- 3. Continuous Jacobian Matrix (F) ---
    % State: [r; v; ba]
    %  d(rDot)/dr = 0,  d(rDot)/dv = I,  d(rDot)/dba = 0
    %  d(vDot)/dr = G,  d(vDot)/dv = 0,  d(vDot)/dba = -DCM_B2ECI
    % d(baDot)/dr = 0, d(baDot)/dv = 0, d(baDot)/dba = 0
    F = [ O3,  I3,        O3;
           G,  O3, -DCM_B2ECI;
          O3,  O3,        O3];
      
    % --- 4. Discrete Approximation ---
    % Phi = I + F*dt (First order Taylor)
    Phi = eye(9) + F * dt;
end