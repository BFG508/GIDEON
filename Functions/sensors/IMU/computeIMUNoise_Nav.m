function [qVRW, qARW] = computeIMUNoise_Nav(IMU)
%==========================================================================
% computeIMUNoise_Nav: Computes raw continuous-time process noise PSD values 
%                      for Accelerometers based on datasheet specs.
%
% Inputs:
%    IMU - IMU structure with fields:
%          .accel.VRW - Velocity Random Walk [μg / sqrt(Hz)]
%          .accel.ARW - Accel Random Walk (Bias Instability) [m/s² / sqrt(hr)]
%
% Outputs:
%    qVRW  - Velocity Random Walk PSD scalar [(m/s²)²/Hz]
%    qARW  - Bias Instability PSD scalar [(m/s²)²/Hz]
%==========================================================================

    g0 = 9.80665; % Standard gravity [m/s²]

    % --- 1. Velocity Random Walk (VRW) ---
    % Physical meaning: White noise on the accelerometer measurement.
    sigmaVRW_metric = (IMU.accel.VRW * 1e-6) * g0; 
    
    % PSD = sigma^2
    qVRW = sigmaVRW_metric^2;
    
    % --- 2. Accel Random Walk (ARW / Bias Instability) ---
    % Physical meaning: Random walk driving the bias drift.
    sigmaARW_metric = IMU.accel.ARW / sqrt(3600);
    
    % PSD = sigma^2
    qARW = sigmaARW_metric^2;
    
end