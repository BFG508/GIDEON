function EKF = updateMAG(EKF, BMeas, B_ECI)
%==========================================================================
% updateMAG: Magnetometer measurement update step for MEKF. Computes Kalman
%            gain from innovation, updates error state and covariance, then
%            applies correction to nominal quaternion and resets error state.
%
% INPUTS:
%    EKF     - EKF structure with fields:
%     .qNom  - Nominal attitude quaternion (4x1), ECI to Body
%     .x     - Error state vector (6x1): [deltaTheta; biasGyro]
%     .P     - Error state covariance matrix (6x6)
%     .R_MAG - Magnetometer measurement noise covariance (3x3) [nT²]
%
%    BMeas   - Measured magnetic field in body frame [nT] (3x1)
%    B_ECI   - Reference magnetic field in ECI frame [nT] (3x1)
%
% OUTPUTS:
%    EKF   - Updated EKF structure
%
% METHOD:
%    1. Predict measurement: BPred = C(qNom) * B_ECI
%    2. Compute innovation: y = BMeas - BPred
%    3. Linearize via Jacobian H = [-[BPred]×, 0]
%    4. Kalman gain: K = P*H' / (H*P*H' + R)
%    5. Update error state: x = x + K*y
%    6. Update covariance using Joseph form for numerical stability
%    7. Apply correction to quaternion and reset attitude error
%==========================================================================

    % --- 1. Predict Measurement ---
    DCM_nom = quat2dcm(EKF.qNom);
    BPred   = DCM_nom * B_ECI;
    
    % --- 2. Innovation (Measurement Residual) ---
    y = BMeas - BPred;
    
    % --- 3. Measurement Jacobian ---
    H = computeMeasJacobian_Att(BPred);
    
    % --- 4. Kalman Gain ---
    S = H * EKF.P * H' + EKF.R_MAG;
    K = EKF.P * H' / S;
    
    % --- 5. Update Error State ---
    EKF.x = EKF.x + K * y;
    
    % --- 6. Update Covariance (Joseph Form) ---
    I     = eye(6);
    IKH   = I - K * H;
    EKF.P = IKH * EKF.P * IKH' + K * EKF.R_MAG * K';
    
    % --- 7. Apply Correction to Nominal Quaternion ---
    deltaQ    = smallAng2quat(EKF.x(1:3));
    EKF.qNom  = quatmultiply(EKF.qNom, deltaQ);
    EKF.qNom  = EKF.qNom / norm(EKF.qNom);
    
    % --- 8. Reset Attitude Error (MEKF Property) ---
    EKF.x(1:3) = zeros(3,1);
end