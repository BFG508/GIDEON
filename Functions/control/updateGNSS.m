function EKF = updateGNSS(EKF, rMeas, vMeas, leverArm, q_ECI2B, omegaBody)
%==========================================================================
% updateGNSS: GNSS PVT measurement update step for the Navigation PV-EKF. 
%             Computes the Kalman gain from the position and velocity 
%             innovation, updates the translational state, and corrects the 
%             accelerometer bias.
%
% INPUTS:
%    EKF       - PV-EKF structure with fields:
%        .x    - State vector (9x1): [r_ECI; v_ECI; biasAccel]
%        .P    - State covariance matrix (9x9)
%        .R    - GNSS measurement noise covariance (6x6)
%
%    rMeas     - Measured antenna position in ECI frame [m] (3x1)
%    vMeas     - Measured antenna velocity in ECI frame [m/s] (3x1)
%    leverArm  - Antenna offset from Center of Mass in Body frame [m] (3x1)
%    q_ECI2B   - Current attitude quaternion (ECI to Body) (4x1)
%    omegaBody - Current angular velocity in Body frame [rad/s] (3x1)
%
% OUTPUTS:
%    EKF       - Updated EKF structure
%
% METHOD:
%    1. Predict antenna measurement using current state and Lever Arm:
%       r_pred = r_CG + DCM_B2I * L_arm
%       v_pred = v_CG + DCM_B2I * (omega x L_arm)
%    2. Compute innovation: y = zMeas - zPred
%    3. Construct linear Jacobian H = [I, 0, 0; 0, I, 0]
%    4. Compute Kalman gain: K = P*H' / (H*P*H' + R)
%    5. Update state vector: x = x + K*y
%    6. Update covariance using Joseph form for numerical stability
%==========================================================================

    % --- 1. Predict Measurement (Including Lever Arm Kinematics) ---
    % Rotation matrix from Body to ECI
    DCM_ECI2B = quat2dcm(q_ECI2B);
    DCM_B2ECI = DCM_ECI2B';
    
    % Antenna position offset in ECI
    r_lever_ECI = DCM_B2ECI * leverArm;
    
    % Antenna velocity offset in ECI (v = w x r)
    v_lever_body = cross(omegaBody, leverArm);
    v_lever_ECI  = DCM_B2ECI * v_lever_body;
    
    % Predicted measurements (Center of Mass + Antenna Offset)
    rPred = EKF.x(1:3) + r_lever_ECI;
    vPred = EKF.x(4:6) + v_lever_ECI;
    
    zPred = [rPred; vPred];
    zMeas = [rMeas; vMeas];
    
    % --- 2. Innovation (Measurement Residual) ---
    y = zMeas - zPred;
    
    % --- 3. Measurement Jacobian (H) ---
    % Since the lever arm depends entirely on attitude (which is an external 
    % input to this translational filter), the partial derivatives of the 
    % measurement w.r.t the state (r, v, bias) are strictly linear.
    % H = d(zPred)/d(x)
    I3 = eye(3);
    O3 = zeros(3,3);
    H  = [ I3, O3, O3;
           O3, I3, O3 ];
       
    % --- 4. Kalman Gain ---
    % S = Innovation Covariance
    S = H * EKF.P * H' + EKF.R_GNSS;
    
    % K = Optimal Gain
    K = EKF.P * H' / S;
    
    % --- 5. Update Error State ---
    EKF.x = EKF.x + K * y;
    
    % --- 6. Update Covariance (Joseph Form) ---
    % Joseph form ensures the covariance matrix remains symmetric and 
    % positive-definite even with numerical round-off errors.
    I9    = eye(9);
    IKH   = I9 - K * H;
    EKF.P = IKH * EKF.P * IKH' + K * EKF.R_GNSS * K';
    
    % Force symmetry explicitly (extra safety net)
    EKF.P = (EKF.P + EKF.P') / 2;
    
end