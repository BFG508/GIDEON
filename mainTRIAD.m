%==========================================================================
% mainTRIAD.m - TRIAD Algorithm for Star Tracker Attitude Determination
%==========================================================================
%
% DESCRIPTION:
%   Complete simulation and validation framework for spacecraft attitude
%   estimation using the QUEST (TRIaxial Attitude Determination) algorithm with
%   multiple star tracker (STR) measurements.
%
% PURPOSE:
%   Demonstrates optimal attitude determination using real Hipparcos stellar
%   catalog and realistic noise models (centroiding, motion blur, mounting
%   misalignment). Validates performance against STR specifications and
%   Cramér-Rao Lower Bound (CRLB).
%
% ALGORITHM:
%   TRIAD solves Wahba's problem using exactly two vector observations by
%   constructing orthonormal triads in both reference (ECI) and body frames.
%   The algorithm:
%     1. Selects primary vector v₁ (highest weight, preserved exactly)
%     2. Selects secondary vector v₂ (second weight, defines perpendicular plane)
%     3. Constructs reference triad: r₁ = v₁, r₂ = v₁×v₂/||v₁×v₂||, r₃ = r₁×r₂
%     4. Constructs body triad: w₁, w₂, w₃ (same procedure)
%     5. Computes DCM directly: A = [w₁ w₂ w₃]·[r₁ r₂ r₃]ᵀ
%  
%   TRIAD is optimal when σ₁ << σ₂ (primary much more accurate than secondary).
%   For equal-weighted observations, QUEST is generally superior.
%
% WORKFLOW:
%   1. Define STR hardware parameters (FOV, resolution, accuracy)
%   2. Load Hipparcos stellar catalog (magnitude ≤ magnitude_limit)
%   3. Set true spacecraft attitude (ground truth)
%   4. Simulate noisy STR measurements
%   5. Execute TRIAD algorithm
%   6. Validate results: errors, CRLB, performance metrics
%   7. Generate diagnostic plots (residuals, quaternion, DCM, vector selection)
%
% KEY FEATURES:
%   - Multi-STR configuration (1 or 2 trackers)
%   - Realistic noise: centroiding, motion blur, misalignment
%   - Real stellar data: Hipparcos catalog (~ 9000 stars)
%   - Comprehensive validation: attitude error, Euler angles, DCM
%   - Performance analysis: CRLB efficiency, dominant noise sources
%   - Publication-quality figures (saved as FIG, PNG, SVG)
%
% OUTPUTS:
%   - qEst:                Optimal attitude quaternion [qw; qx; qy; qz]
%   - DCMEst:              Direction Cosine Matrix (3×3)
%   - Attitude errors:     Total, Euler angles, DCM Frobenius norm
%   - Performance metrics: CRLB, efficiency, residual statistics
%   - Diagnostic figures:  5 plots in Figures/ directory
%
% CONFIGURATION:
%   Edit Section 1 (initializeSTR.m) for:
%   - STR hardware specs
%   - Mounting misalignment
%   - Number of STRs
%   - Noise levels
%
%   Edit Section 3 for:
%   - Spacecraft true attitude
%   - Spacecraft angular velocity

% TRIAD vs. QUEST:
% TRIAD advantages:
%   - Simpler computation (no eigenvalue problem)
%   - Optimal when primary vector is much more accurate
%   - Deterministic solution (no iterative convergence)
%
% TRIAD limitations:
%   - Uses only 2 vectors (discards additional observations)
%   - Asymmetric treatment (primary preserved, secondary projected)
%   - Suboptimal for equal-weighted observations
%
% Recommendation: Use TRIAD when σ₁ << σ₂ (weight ratio > 3:1).
%                 Otherwise, use QUEST for better multi-star utilization.
%==========================================================================

clear;
close all;
clc;

%% ========================================================================
%  1. STAR TRACKER PARAMETERS (typical values from commercial datasheets)
% =========================================================================

nSTR = 2;
STR  = initializeSTR(nSTR);

%% ========================================================================
%  2. LOAD AND FILTER HIPPARCOS STELLAR CATALOG
% =========================================================================

% Load real Hipparcos catalog (cached in Data/ folder for speed)
catalog = loadHipparcosStellarCatalog(STR(1).magnitudeLimit);

%% ========================================================================
%  3. DEFINE TRUE ATTITUDE (GROUND TRUTH)
% =========================================================================

% True attitude quaternion (scalar-first convention: [qs; qv])
angleTrue = rand * 2*pi;                   % Rotation angle [rad]
 axisTrue = randn(3,1);
 axisTrue = axisTrue / norm(axisTrue);     % Unit rotation axis
 
    qTrue = [cos(angleTrue/2); ...         % Scalar part
             sin(angleTrue/2) * axisTrue]; % Vector part

% Convert quaternion to DCM (ECI-to-Body transformation)
DCMTrue = quat2dcm(qTrue);

% Satellite angular velocity in body frame [rad/s]
omegaBody = randn(3,1) * 1e-5;

%% ========================================================================
%  4. GENERATE SYNTHETIC STAR TRACKER MEASUREMENTS
% =========================================================================

% Simulate observations from all configured star trackers
meas = generateSTRMeasurements(STR, nSTR, catalog, DCMTrue, omegaBody);

%% ========================================================================
%  5. TRIAD ALGORITHM - ATTITUDE ESTIMATION
% =========================================================================

% Solve for attitude using TRIAD algorithm
[qEst, DCMEst, triadInfo] = solveTRIADAttitude(meas, nSTR);

%% ========================================================================
% 6. VALIDATION AND ERROR ANALYSIS
% =========================================================================

fprintf('\n=== Attitude Error Analysis ===\n');

% Quaternion Error
qErr     = quatmultiply(qEst, quatinv(qTrue));
angleErr = rad2deg(2 * acos(min(abs(qErr(1)), 1))) * 3600; % Scalar-first: q(1)

fprintf('\nQuaternion Comparison:\n');
fprintf('  qTrue: [%.6f, %.6f, %.6f, %.6f]\n', qTrue);
fprintf('  qEst:  [%.6f, %.6f, %.6f, %.6f]\n', qEst);
fprintf('  qErr:  [%.6f, %.6f, %.6f, %.6f]\n', qErr);

fprintf('\nTotal Angular Error:\n');
fprintf('  Magnitude:   %.4f deg = %.2f arcsec  ', angleErr, angleErr);
if angleErr < 5
    fprintf('✓ Excellent\n');
elseif angleErr < 20
    fprintf('✓ Good\n');
elseif angleErr < 100
    fprintf('○ Acceptable\n');
else
    fprintf('⚠ High error\n');
end

% Euler Angles Error (ZYX: Yaw-Pitch-Roll)
eulerTrue = rad2deg(dcm2euler_ZYX(DCMTrue));
eulerEst  = rad2deg(dcm2euler_ZYX(DCMEst));
eulerErr  = eulerEst - eulerTrue;

fprintf('\nEuler Angles (ZYX: Yaw-Pitch-Roll):\n');
fprintf('  Component    True        Estimated   Error       Error\n');
fprintf('               [deg]         [deg]     [deg]      [arcsec]\n');
fprintf('  ─────────────────────────────────────────────────────────\n');
fprintf('  Yaw   (Z)    %7.3f     %7.3f    %+7.4f   %+8.2f\n', ...
        eulerTrue(1), eulerEst(1), eulerErr(1), eulerErr(1)*3600);
fprintf('  Pitch (Y)    %7.3f     %7.3f    %+7.4f   %+8.2f\n', ...
        eulerTrue(2), eulerEst(2), eulerErr(2), eulerErr(2)*3600);
fprintf('  Roll  (X)    %7.3f     %7.3f    %+7.4f   %+8.2f\n', ...
        eulerTrue(3), eulerEst(3), eulerErr(3), eulerErr(3)*3600);

[maxEulerErr, maxAxis] = max(abs(eulerErr));
axisNames = {'Yaw', 'Pitch', 'Roll'};
fprintf('\n  Max Error:   %s Axis (%.2f arcsec)\n', ...
        axisNames{maxAxis}, maxEulerErr * 3600);

% Comparison with STR Specifications
fprintf('\n--- Performance vs. STR Specifications ---\n');

    % Cross-boresight accuracy (affects Yaw/Pitch)
    crossboreSpec_3sigma = 3 * STR(1).attitudeAccuracyCrossbore;
    crossboreERR         = sqrt(eulerErr(1)^2 + eulerErr(2)^2) * 3600;
    
    fprintf('  Cross-boresight:\n');
    fprintf('    Spec (1σ):       %6.1f arcsec\n', STR(1).attitudeAccuracyCrossbore);
    fprintf('    Spec (3σ):       %6.1f arcsec\n', crossboreSpec_3sigma);
    fprintf('    Measured error:  %6.2f arcsec  ', crossboreERR);
    if crossboreERR < STR(1).attitudeAccuracyCrossbore
        fprintf('✓ Better than 1σ\n');
    elseif crossboreERR < crossboreSpec_3sigma
        fprintf('✓ Within 3σ\n');
    else
        fprintf('⚠ Exceeds 3σ spec\n');
    end
    
    % Roll accuracy (affects Roll axis)
    rollSpec_3sigma = 3 * STR(1).attitudeAccuracyRoll;
    rollErr         = abs(eulerErr(3)) * 3600;
    
    fprintf('  Roll-axis:\n');
    fprintf('    Spec (1σ):       %6.1f arcsec\n', STR(1).attitudeAccuracyRoll);
    fprintf('    Spec (3σ):       %6.1f arcsec\n', rollSpec_3sigma);
    fprintf('    Measured error:  %6.2f arcsec  ', rollErr);
    if rollErr < STR(1).attitudeAccuracyRoll
        fprintf('✓ Better than 1σ\n');
    elseif rollErr < rollSpec_3sigma
        fprintf('✓ Within 3σ\n');
    else
        fprintf('⚠ Exceeds 3σ spec\n');
    end
    
    % Overall attitude error vs. cross-boresight spec
    fprintf('  Overall attitude:\n');
    fprintf('    Measured error:  %6.2f arcsec  ', angleErr);
    if angleErr < STR(1).attitudeAccuracyCrossbore
        fprintf('✓ Better than 1σ spec\n');
    elseif angleErr < crossboreSpec_3sigma
        fprintf('✓ Within 3σ spec\n');
    else
        fprintf('⚠ Exceeds 3σ spec\n');
    end
    
% Cramér-Rao Lower Bound (Theoretical Best Performance)
fprintf('\n--- Theoretical Performance Limits (CRLB) ---\n');
    % Single-axis pointing precision from centroiding
    sigmaCentroid = rad2deg(  STR(1).centroidAccuracy * STR(1).pixelSize ... 
                            / STR(1).focalLength) * 3600;
    
    fprintf('  Centroiding noise:   %.3f pixels (1σ)\n', STR(1).centroidAccuracy);
    fprintf('  Single-star CRLB:    %.2f arcsec\n'     , sigmaCentroid);
    
    % Multi-star improvement: σ_combined ≈ σ_single / sqrt(N)
    nEffective = triadInfo.nObs;
    sigmaCRLB  = sigmaCentroid / sqrt(nEffective);
    
    fprintf('  Multi-star CRLB:     %.2f arcsec (N = %d stars)\n', ...
            sigmaCRLB, nEffective);
    
    % Compare actual error with theoretical limit
    fprintf('  Measured error:      %.2f arcsec\n', angleErr);
    fprintf('  Ratio (error/CRLB):  %.2fx  '      , angleErr / sigmaCRLB);
    if angleErr < 2 * sigmaCRLB
        fprintf('✓ Within 2σ of CRLB\n');
    elseif angleErr < 3 * sigmaCRLB
        fprintf('○ Within 3σ of CRLB\n');
    else
        fprintf('⚠ Exceeds 3σ (systematic error present)\n');
    end
    
    % Identify dominant noise source
    if angleErr < 2 * sigmaCRLB
        dominantNoise = 'Centroid noise (measurement-limited)';
    elseif triadInfo.residualMean > 100
        dominantNoise = 'Motion blur (dynamics-limited)';
    elseif triadInfo.residualStd / triadInfo.residualMean > 0.8
        dominantNoise = 'Outliers or inconsistent measurements';
    else
        dominantNoise = 'Mixed noise sources';
    end
    fprintf('  Dominant source:     %s\n', dominantNoise);

% DCM Error Matrix
fprintf('\n--- Direction Cosine Matrix Error ---\n');

DCMErr = DCMEst * DCMTrue'; % Should be close to identity
DCM_errAngles = rad2deg([asin(DCMErr(3,2) - DCMErr(2,3)) / 2;
                         asin(DCMErr(1,3) - DCMErr(3,1)) / 2;
                         asin(DCMErr(2,1) - DCMErr(1,2)) / 2]) * 3600;

fprintf('  Small-angle errors:  [%.2f, %.2f, %.2f] arcsec\n', DCM_errAngles);
fprintf('  Frobenius norm:      %.2e  ', norm(DCMErr - eye(3), 'fro'));
if norm(DCMErr - eye(3), 'fro') < 1e-5
    fprintf('✓ Excellent\n');
elseif norm(DCMErr - eye(3), 'fro') < 1e-3
    fprintf('✓ Good\n');
else
    fprintf('○ Acceptable\n');
end

% TRIAD-Specific Performance Analysis
fprintf('\n--- TRIAD-Specific Analysis ---\n');
    % TRIAD uses exactly 2 vectors (primary and secondary)
    % Theoretical accuracy depends on vector weights and separation angle
    nVectorsUsed  = 2;
    weightRatio   = triadInfo.weightPrimary / triadInfo.weightSecondary;
    separationAng = deg2rad(triadInfo.angleBetween);
    
    % TRIAD CRLB approximation (from Shuster 2007):
    % σ_TRIAD² ≈ (σ₁² + σ₂²·sin²θ) / sin²θ
    % where θ is angle between vectors, σ₁ and σ₂ are measurement uncertainties
    % Simplified: σ_TRIAD ≈ σ_centroid / (sin(θ) * sqrt(2))
    sigmaTriadFactor = 1 / (sin(separationAng) * sqrt(2));
    sigmaTriad       = sigmaCentroid * sigmaTriadFactor;
    
    fprintf('  TRIAD CRLB:          %.2f arcsec (2 vectors, θ = %.1f°)\n', ...
            sigmaTriad, triadInfo.angleBetween);
    fprintf('  Primary weight:      %.4f (most accurate vector)\n', triadInfo.weightPrimary);
    fprintf('  Secondary weight:    %.4f (ratio = %.2f:1)\n', ...
            triadInfo.weightSecondary, weightRatio);
    
    % Compare actual error with TRIAD theoretical limit
    fprintf('  Measured error:      %.2f arcsec\n', angleErr);
    fprintf('  Ratio (error/TRIAD): %.2fx  ', angleErr / sigmaTriad);
    if angleErr < 2 * sigmaTriad
        fprintf('✓ Within 2σ of TRIAD CRLB\n');
    elseif angleErr < 3 * sigmaTriad
        fprintf('○ Within 3σ of TRIAD CRLB\n');
    else
        fprintf('⚠ Exceeds 3σ (systematic error present)\n');
    end

% TRIAD Optimality Assessment
fprintf('\n--- TRIAD Optimality Assessment ---\n');
    % Weight distribution analysis
    if weightRatio > 3.0
        fprintf('  Weight distribution: ✓ TRIAD optimal (σ₁ << σ₂, ratio %.1f:1)\n', weightRatio);
        fprintf('                       Primary vector dominates accuracy as expected.\n');
    elseif weightRatio > 1.5
        fprintf('  Weight distribution: ○ TRIAD acceptable (moderate ratio %.1f:1)\n', weightRatio);
        fprintf('                       QUEST may provide similar or slightly better accuracy.\n');
    else
        fprintf('  Weight distribution: ⚠ Vectors have similar weights (%.1f:1)\n', weightRatio);
        fprintf('                       QUEST recommended: Better utilizes all %d observations.\n', triadInfo.nObs);
    end
    
    % Vector separation analysis
    if triadInfo.angleBetween > 30 && triadInfo.angleBetween < 150
        fprintf('  Vector separation:   ✓ Good conditioning (%.1f°)\n', triadInfo.angleBetween);
    elseif triadInfo.angleBetween > 15 && triadInfo.angleBetween < 165
        fprintf('  Vector separation:   ○ Acceptable (%.1f°)\n', triadInfo.angleBetween);
    else
        fprintf('  Vector separation:   ⚠ Poor conditioning (%.1f°, near-collinear)\n', ...
                triadInfo.angleBetween);
    end
    
    % Unused observations warning
    if triadInfo.nObs > 2
        fprintf('  Additional stars:    %d observations not used by TRIAD\n', ...
                triadInfo.nObs - 2);
        fprintf('                       Consider QUEST to utilize all available measurements.\n');
    end

%% ========================================================================
%  7. RESULTS VISUALIZATION
% =========================================================================

% Generate comprehensive validation plots
saveFlag = 1;
plotTRIADResults(meas, nSTR, STR, DCMTrue, DCMEst, ...
                 qTrue, qEst, triadInfo, saveFlag);

fprintf('\n=== All tasks completed successfully ===\n');