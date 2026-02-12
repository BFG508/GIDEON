function R_MAG = computeMAGNoise(MAG)
%==========================================================================
% computeMAGNoise: Computes MAG measurement noise covariance matrix based 
%                  on white noise characteristics. Assumes isotropic,
%                  uncorrelated Gaussian noise across sensor axes.
%
% Inputs:
%    MAG   - Magnetometer structure with fields:
%            .sigmaWhite - White noise RMS [nT] (1σ)
%
% Outputs:
%    R_MAG - Measurement noise covariance matrix (3x3) [nT²]
%            Diagonal matrix representing isotropic field measurement noise
%
% Method:
%    White noise from sensor datasheet (noise density × √bandwidth)
%    Bias instability is modeled in process noise, not measurement noise
%    Quantization noise is negligible for modern ADCs (>12 bits)
%==========================================================================

    % --- 1. White Noise RMS ---
    sigmaWhite = MAG.sigmaWhite;
    
    % --- 2. Isotropic Measurement Noise Covariance ---
    R_MAG = sigmaWhite^2 * eye(3);

end
