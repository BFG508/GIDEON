function R_STR = computeSTRNoise(STR)
%==========================================================================
% computeSTRNoise: Computes Star Tracker measurement noise covariance matrix
%                  based on centroiding error, multi-star averaging, and
%                  boresight geometry.
%
% Inputs:
%    STR - Star Tracker structure array with fields:
%        .centroidAccuracy - Centroiding error [pixels] (1σ)
%        .pixelSize        - Pixel pitch            [m]
%        .focalLength      - Focal length           [m]
%        .FOV              - Field of view        [deg]
%        .magnitudeLimit   - Magnitude limit
%
% Outputs:
%    R_STR - Measurement noise covariance matrix (3x3) [rad²]
%            Represents attitude error uncertainty in body frame
%            Anisotropic: assumes STR boresight along +Z body axis
%
% Method:
%    1. Convert pixel centroiding error to angular error per star
%    2. Estimate number of tracked stars from FOV and magnitude limit
%    3. Apply √N averaging improvement for multi-star solution
%    4. Account for boresight vs roll-axis accuracy difference
%    5. Build anisotropic covariance (cross-boresight better than roll)
%
% Reference:
%    Liebe, C. C. (2002). 
%    "Accuracy Performance of Star Trackers - A Tutorial".
%    IEEE Transactions on Aerospace and Electronic Systems, 38(2), 587-599.
%==========================================================================

    % --- 1. Single-Star Angular Noise ---
    sigmaCentroidPix = STR(1).centroidAccuracy;
    sigmaStarRad     = sigmaCentroidPix * STR(1).pixelSize / STR(1).focalLength;
    
    % --- 2. Estimate Number of Tracked Stars ---
    if isfield(STR, 'FOV') && isfield(STR, 'magnitudeLimit')
        stellarDensity = 600 * 10^(0.6 * (STR(1).magnitudeLimit - 6.0));
        fovSteradian   = (pi * (STR(1).FOV/2)^2) / (180^2);
        nStarsAvg      = max(5, stellarDensity * fovSteradian * 0.7);
    else
        nStarsAvg = 10;
    end
    
    % --- 3. Multi-Star Averaging ---
    sigmaCrossbore = sigmaStarRad / sqrt(nStarsAvg);
    
    % --- 4. Roll-Axis Degradation ---
    sigmaRoll = 4 * sigmaCrossbore;
    
    % --- 5. Anisotropic Measurement Noise Covariance ---
    R_STR = diag([sigmaCrossbore^2, sigmaCrossbore^2, sigmaRoll^2]);
end
