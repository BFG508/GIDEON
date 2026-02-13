function EKF = predictEKF_Att(EKF, omegaMeas, dt)
%==========================================================================
% predictEKF_Att: Propagates MEKF state and covariance through the prediction
%                 step using gyrp measurements. Updates the nominal quaternion, 
%                 propagates the error-state covariance, and resets the error state.
%
% Inputs:
%    EKF       - EKF structure with fields:
%        .qNom - Nominal attitude quaternion (4x1), ECI to Body frame
%        .x    - Error state vector (6x1): [deltaTheta; biasGyro]
%        .P    - Error state covariance matrix (6x6)
%        .Q    - Process noise covariance matrix (6x6)
%
%    omegaMeas - Measured angular velocity from gyroscope [rad/s] (3x1)
%    dt        - Time step [s]
%
% Outputs:
%    ekf        - Updated EKF structure propagated.
%==========================================================================

    % --- 1. Propagate nominal quaternion with bias-corrected gyro ---
    omegaCorrected = omegaMeas - EKF.x(4:6);
    EKF.qNom       = quatint(EKF.qNom, omegaCorrected, dt);
    
    % --- 2. Propagate error covariance: P = Φ*P*Φ' + Q ---
    Phi   = computeSTM(omegaCorrected, dt);
    EKF.P = Phi * EKF.P * Phi' + EKF.Q * dt;
    EKF.P = (EKF.P + EKF.P') / 2;
end
