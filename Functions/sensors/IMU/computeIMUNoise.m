function Q = computeIMUNoise(IMU)
%==========================================================================
% computeIMUNoise: Computes discrete-time process noise covariance matrix
%                  for MEKF error state propagation. Models gyroscope
%                  measurement noise and bias random walk.
%
% Inputs:
%    IMU - IMU structure with fields:
%          .gyro.ARW - Angle Random Walk [deg/√h]
%          .gyro.RRW - Rate Random Walk [deg/h/√h]
%
% Outputs:
%    Q   - Process noise spectral density matrix (6x6) [(rad)²/s]
%          Block-diagonal structure: Q = blkdiag(QAtt, QBias)
%          where QAtt: attitude error noise (3x3)
%                QBias: gyro bias drift noise (3x3)
%
% Method:
%    Continuous-time spectral density: Qc = blkdiag(σ_ARW²·I₃, σ_RRW²·I₃)
%    Discrete-time approximation: Qd ≈ Qc · dt
%
% REFERENCE:
%    IEEE Standard Specification Format Guide and Test Procedure for Linear (1998).
%==========================================================================

    % --- 1. Attitude Error Process Noise (Angle Random Walk) ---
    sigmaARW = (deg2rad(IMU.gyro.ARW) / sqrt(3600))^2;
    QAtt     = sigmaARW * eye(3);
    
    % --- 2. Gyro Bias Process Noise (Rate Random Walk) ---
    sigmaRRW = (deg2rad(IMU.gyro.RRW) / 3600 / sqrt(3600))^2;
    QBias    = sigmaRRW * eye(3);
    
    % --- 3. Block-Diagonal Process Noise Matrix ---
    Q        = blkdiag(QAtt, QBias);
    
end
