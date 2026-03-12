function EKF = updateGNSS(EKF, rMeas, vMeas, leverArm, q_ECI2B, omegaBody)
%==========================================================================
% updateGNSS: GNSS PVT measurement update step for the Navigation PV-EKF. 
%             Computes the Kalman gain from the position and velocity 
%             innovation, updates the translational state, and corrects the 
%             accelerometer bias.
%
% Inputs:
%    EKF       - PV-EKF structure with fields:
%        .x    - State vector: [r_ECI; v_ECI; biasAccel; rGM; vGM]     , (15x1)
%        .P    - State covariance matrix                               , (15x15)
%        .R    - GNSS measurement noise covariance                     , (6x6)
%    rMeas     - Measured antenna position in ECI frame             [m], (3x1)
%    vMeas     - Measured antenna velocity in ECI frame           [m/s], (3x1)
%    leverArm  - Antenna offset from Center of Mass in Body frame   [m], (3x1)
%    q_ECI2B   - Current attitude quaternion (ECI to Body)             , (4x1)
%    omegaBody - Current angular velocity in Body frame         [rad/s], (3x1)
%
% Outputs:
%    EKF       - Updated EKF structure
%
% Method:
%    1. Predict antenna measurement using current state, Lever Arm, and GM:
%       rPred = rCG + DCM_B2I * LArm + rGM
%       vPred = vCG + DCM_B2I * (omega x LArm) + vGM
%    2. Compute innovation: y = zMeas - zPred
%    3. Construct linear Jacobian H = [I, 0, 0, I, 0; 0, I, 0, 0, I]
%    4. Compute Kalman gain: K = P*H' / (H*P*H' + R)
%    5. Update state vector: x = x + K*y
%    6. Update covariance using Joseph form for numerical stability
%==========================================================================

    % --- 1. Predict Measurement ---
        % Rotation matrix from Body to ECI
        DCM_ECI2B = quat2dcm(q_ECI2B);
        DCM_B2ECI = DCM_ECI2B';
        
        % Antenna position offset in ECI
        rLever_ECI = DCM_B2ECI * leverArm;
        
        % Antenna velocity offset in ECI (v = w x r)
        vLever_body = cross(omegaBody, leverArm);
        vLever_ECI  = DCM_B2ECI * vLever_body;
        
        % Predicted measurements (Center of Mass + Antenna Offset + GM Error)
        rPred = EKF.x(1:3) + rLever_ECI + EKF.x(10:12);
        vPred = EKF.x(4:6) + vLever_ECI + EKF.x(13:15);
        
        zPred = [rPred; vPred];
        zMeas = [rMeas; vMeas];
    
    % --- 2. Innovation (Measurement Residual) ---
    y = zMeas - zPred;
    
    % --- 3. Measurement Jacobian (H) ---
    % Since the lever arm depends entirely on attitude (which is an external 
    % input to this translational filter), the partial derivatives of the 
    % measurement w.r.t the state (r, v, bias, rGM, vGM) are strictly linear.
    % H = d(zPred)/d(x)
    I3 = eye(3);
    O3 = zeros(3,3);
    H  = [I3, O3, O3, I3, O3;
          O3, I3, O3, O3, I3];
       
    % --- 4. Kalman Gain ---
    S = H * EKF.P * H' + EKF.R_GNSS;
    K = EKF.P * H' / S;
    
    % --- 5. Update Error State ---
    EKF.x = EKF.x + K * y;
    
    % --- 6. Update Covariance (Joseph Form) ---
    I15   = eye(15);
    IKH   = I15 - K * H;
    EKF.P = IKH * EKF.P * IKH' + K * EKF.R_GNSS * K';
    EKF.P = (EKF.P + EKF.P') / 2;
    
end