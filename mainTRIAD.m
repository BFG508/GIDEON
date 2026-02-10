%==========================================================================
% mainTRIAD.m - TRIAD Algorithm for Star Tracker Attitude Determination
%==========================================================================
%
% DESCRIPTION:
%   Complete simulation and validation framework for spacecraft attitude
%   estimation using the TRIAD (TRIaxial Attitude Determination) algorithm
%   with multiple star tracker (STR) measurements.
%
% PURPOSE:
%   Demonstrates two-vector attitude determination using real Hipparcos
%   stellar catalog and realistic noise models (centroiding, motion blur,
%   mounting misalignment). Validates performance against STR specifications
%   and theoretical Cramér-Rao Lower Bound (CRLB) for TRIAD geometry.
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
%   2. Load Hipparcos stellar catalog (magnitude ≤ 6.5)
%   3. Set true spacecraft attitude (ground truth)
%   4. Simulate noisy STR measurements (up to 2 STRs)
%   5. Execute TRIAD algorithm (two-vector deterministic method)
%   6. Validate results: errors, CRLB, performance metrics
%   7. Generate diagnostic plots (residuals, quaternion, DCM, vector selection)
%
% KEY FEATURES:
%   - Multi-STR configuration (1 or 2 trackers)
%   - Realistic noise: centroiding (0.15 px), motion blur, misalignment
%   - Real stellar data: Hipparcos catalog (~9000 stars)
%   - Comprehensive validation: attitude error, Euler angles, DCM
%   - TRIAD optimality assessment: weight ratio, angular separation
%   - Performance analysis: CRLB efficiency, dominant noise sources
%   - Publication-quality figures (saved as FIG, PNG, SVG)
%
% OUTPUTS:
%   - q_estimated: Optimal attitude quaternion [qw; qx; qy; qz]
%   - DCM_estimated: Direction Cosine Matrix (3×3)
%   - Attitude errors: Total, Euler angles, DCM Frobenius norm
%   - Performance metrics: CRLB, efficiency, residual statistics
%   - Diagnostic figures: 6 plots in Figures/ directory
%
% CONFIGURATION:
%   Edit Section 1 for:
%   - STR hardware specs (FOV, resolution, focal length)
%   - Mounting misalignment (roll, pitch, yaw errors)
%   - Number of STRs (N_STR = 1 or 2)
%   - Noise levels (centroiding accuracy, angular rate)
%
%   Edit Section 3 for:
%   - True attitude (quaternion or Euler angles)
%   - Spacecraft angular velocity (omega_body)

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

% Mounting misalignment errors (typical mechanical tolerances)
roll_error  =  0.000; % Mounting error around X-axis [deg]
pitch_error = -0.000; % Mounting error around Y-axis [deg]
yaw_error   =  0.000; % Mounting error around Z-axis [deg]

% Construct misalignment DCM (small perturbation from ideal alignment)
DCM_misalignment = angle2dcm(deg2rad(yaw_error), ...
                             deg2rad(pitch_error), ...
                             deg2rad(roll_error), ...
                             'ZYX');

% Base STR specifications (representative of high-accuracy units)
STR.FOV_deg                     = 20;           % Field of view               [deg]
STR.resolution                  = [2048, 2048]; % Sensor resolution           [pixels]
STR.pixel_size                  = 5.5e-6;       % Pixel pitch                 [m]
STR.focal_length                = 0.065;        % Optical focal length        [m]
STR.centroid_accuracy           = 0.15;         % Centroiding precision       [pixels, 1σ]
STR.magnitude_limit             = 6.5;          % Visual magnitude threshold  [-]
STR.update_rate                 = 4;            % Measurement frequency       [Hz]
STR.attitude_accuracy_crossbore = 5;            % Cross-boresight accuracy    [arcsec, 1σ]
STR.attitude_accuracy_roll      = 20;           % Roll-axis accuracy          [arcsec, 1σ]
STR.max_angular_rate            = 1.0;          % Maximum trackable rate      [deg/s]

% Number of star trackers in configuration
N_STR = 2; % 1 or 2 STRs supported

% STR1: Nominal boresight aligned with +Z body axis
STR(1).DCM_B2S                     = DCM_misalignment;
STR(1).FOV_deg                     = STR.FOV_deg;
STR(1).resolution                  = STR.resolution;
STR(1).pixel_size                  = STR.pixel_size;
STR(1).focal_length                = STR.focal_length;
STR(1).centroid_accuracy           = STR.centroid_accuracy;
STR(1).magnitude_limit             = STR.magnitude_limit;
STR(1).update_rate                 = STR.update_rate;
STR(1).attitude_accuracy_crossbore = STR.attitude_accuracy_crossbore;
STR(1).attitude_accuracy_roll      = STR.attitude_accuracy_roll;

% STR2: 90° rotation around Y-axis for lateral coverage (if N_STR = 2)
if N_STR == 2
    % Compose misalignment with nominal 90° Y-rotation
    DCM_nominal_90Y = [ cos(pi/2), 0, sin(pi/2); 
                                0, 1,         0; 
                       -sin(pi/2), 0, cos(pi/2)];
    
    STR(2).DCM_B2S                     = DCM_misalignment * DCM_nominal_90Y;
    STR(2).FOV_deg                     = STR.FOV_deg;
    STR(2).resolution                  = STR.resolution;
    STR(2).pixel_size                  = STR.pixel_size;
    STR(2).focal_length                = STR.focal_length;
    STR(2).centroid_accuracy           = STR.centroid_accuracy;
    STR(2).magnitude_limit             = STR.magnitude_limit;
    STR(2).update_rate                 = STR.update_rate;
    STR(2).attitude_accuracy_crossbore = STR.attitude_accuracy_crossbore;
    STR(2).attitude_accuracy_roll      = STR.attitude_accuracy_roll;
end

%% ========================================================================
%  2. LOAD AND FILTER HIPPARCOS STELLAR CATALOG
% =========================================================================

% Load real Hipparcos catalog (cached in Data/ folder for speed)
catalog = loadHipparcosStellarCatalog(STR(1).magnitude_limit);

%% ========================================================================
%  3. DEFINE TRUE ATTITUDE (GROUND TRUTH)
% =========================================================================

% True attitude quaternion (scalar-first convention: [qs; qv])
% Rotation: 30° about normalized axis [1,1,1]
angle_true = deg2rad(30);                 % Rotation angle [rad]
axis_true = [1; 1; 1] / sqrt(3);          % Unit rotation axis
q_true = [cos(angle_true/2); ...          % Scalar part
          sin(angle_true/2) * axis_true]; % Vector part

% Convert quaternion to DCM (ECI-to-Body transformation)
DCM_true = quat2dcm_custom(q_true);

% Satellite angular velocity in body frame [rad/s]
omega_body = [ 0.000; 
               0.000; 
              -0.000];

%% ========================================================================
%  4. GENERATE SYNTHETIC STAR TRACKER MEASUREMENTS
% =========================================================================

% Simulate observations from all configured star trackers
meas = generateSTRMeasurements(STR, N_STR, catalog, DCM_true, omega_body);

%% ========================================================================
%  5. TRIAD ALGORITHM - ATTITUDE ESTIMATION
% =========================================================================

% Solve for attitude using TRIAD algorithm
[q_estimated, DCM_estimated, triad_info] = solveTRIADAttitude(meas, N_STR);

%% ========================================================================
% 6. VALIDATION AND ERROR ANALYSIS
% =========================================================================
fprintf('\n=== Attitude Error Analysis ===\n');

% Quaternion Error
q_error = quatmultiply_custom(q_estimated, quatinv_custom(q_true));
angle_error_rad = 2 * acos(min(abs(q_error(1)), 1)); % Scalar-first: q(1)
angle_error_deg = rad2deg(angle_error_rad);
angle_error_arcsec = angle_error_deg * 3600;

fprintf('\nQuaternion Comparison:\n');
fprintf(' q_true:      [%.6f, %.6f, %.6f, %.6f]\n', q_true);
fprintf(' q_estimated: [%.6f, %.6f, %.6f, %.6f]\n', q_estimated);
fprintf(' q_error:     [%.6f, %.6f, %.6f, %.6f]\n', q_error);

fprintf('\nTotal Angular Error:\n');
fprintf(' Magnitude: %.4f deg = %.2f arcsec ', angle_error_deg, angle_error_arcsec);
if angle_error_arcsec < 5
    fprintf('✓ Excellent\n');
elseif angle_error_arcsec < 20
    fprintf('✓ Good\n');
elseif angle_error_arcsec < 100
    fprintf('○ Acceptable\n');
else
    fprintf('⚠ High error\n');
end

% Euler Angles Error (ZYX: Yaw-Pitch-Roll)
euler_true = rad2deg(dcm2euler_ZYX(DCM_true));
euler_estimated = rad2deg(dcm2euler_ZYX(DCM_estimated));
euler_error = euler_estimated - euler_true;

fprintf('\nEuler Angles (ZYX: Yaw-Pitch-Roll):\n');
fprintf('  Component    True      Estimated   Error      Error\n');
fprintf('               [deg]       [deg]     [deg]     [arcsec]\n');
fprintf(' ─────────────────────────────────────────────────────────\n');
fprintf('  Yaw   (Z)  %7.3f     %7.3f    %+7.4f  %+8.2f\n', ...
    euler_true(1), euler_estimated(1), euler_error(1), euler_error(1)*3600);
fprintf('  Pitch (Y)  %7.3f     %7.3f    %+7.4f  %+8.2f\n', ...
    euler_true(2), euler_estimated(2), euler_error(2), euler_error(2)*3600);
fprintf('  Roll  (X)  %7.3f     %7.3f    %+7.4f  %+8.2f\n', ...
    euler_true(3), euler_estimated(3), euler_error(3), euler_error(3)*3600);

[max_euler_error, max_axis] = max(abs(euler_error));
axis_names = {'Yaw', 'Pitch', 'Roll'};
fprintf('\n  Max error: %s axis (%.2f arcsec)\n', ...
    axis_names{max_axis}, max_euler_error * 3600);

% Comparison with STR Specifications
fprintf('\n--- Performance vs. STR Specifications ---\n');

    % Cross-boresight accuracy (affects Yaw/Pitch)
    crossbore_spec_3sigma = 3 * STR(1).attitude_accuracy_crossbore;
    crossbore_error = sqrt(euler_error(1)^2 + euler_error(2)^2) * 3600;
    fprintf('  Cross-boresight:\n');
    fprintf('    Spec (1σ):      %6.1f arcsec\n', STR(1).attitude_accuracy_crossbore);
    fprintf('    Spec (3σ):      %6.1f arcsec\n', crossbore_spec_3sigma);
    fprintf('    Measured error: %6.2f arcsec ', crossbore_error);
    if crossbore_error < STR(1).attitude_accuracy_crossbore
        fprintf('✓ Better than 1σ\n');
    elseif crossbore_error < crossbore_spec_3sigma
        fprintf('✓ Within 3σ\n');
    else
        fprintf('⚠ Exceeds 3σ spec\n');
    end
    
    % Roll accuracy (affects Roll axis)
    roll_spec_3sigma = 3 * STR(1).attitude_accuracy_roll;
    roll_error = abs(euler_error(3)) * 3600;
    fprintf('  Roll-axis:\n');
    fprintf('    Spec (1σ):      %6.1f arcsec\n', STR(1).attitude_accuracy_roll);
    fprintf('    Spec (3σ):      %6.1f arcsec\n', roll_spec_3sigma);
    fprintf('    Measured error: %6.2f arcsec ', roll_error);
    if roll_error < STR(1).attitude_accuracy_roll
        fprintf('✓ Better than 1σ\n');
    elseif roll_error < roll_spec_3sigma
        fprintf('✓ Within 3σ\n');
    else
        fprintf('⚠ Exceeds 3σ spec\n');
    end
    
    % Overall attitude error vs. cross-boresight spec
    fprintf('  Overall attitude:\n');
    fprintf('    Measured error: %6.2f arcsec ', angle_error_arcsec);
    if angle_error_arcsec < STR(1).attitude_accuracy_crossbore
        fprintf('✓ Better than 1σ spec\n');
    elseif angle_error_arcsec < crossbore_spec_3sigma
        fprintf('✓ Within 3σ spec\n');
    else
        fprintf('⚠ Exceeds 3σ spec\n');
    end
    
    % Cramér-Rao Lower Bound (Theoretical Best Performance for TRIAD)
    fprintf('\n--- Theoretical Performance Limits (CRLB for TRIAD) ---\n');
    
    % Single-axis pointing precision from centroiding
    sigma_centroid_rad = STR(1).centroid_accuracy * STR(1).pixel_size / STR(1).focal_length;
    sigma_centroid_arcsec = rad2deg(sigma_centroid_rad) * 3600;
    fprintf('  Centroiding noise: %.3f pixels (1σ)\n', STR(1).centroid_accuracy);
    fprintf('  Single-star CRLB:  %.2f arcsec\n', sigma_centroid_arcsec);

% TRIAD uses exactly 2 vectors (primary and secondary)
% Theoretical accuracy depends on vector weights and separation angle
N_vectors_used = 2;
weight_ratio = triad_info.vector1_weight / triad_info.vector2_weight;
angle_sep_rad = deg2rad(triad_info.angle_between_vectors);

% TRIAD CRLB approximation (from Shuster 2007):
% σ_TRIAD² ≈ (σ₁² + σ₂²·sin²θ) / sin²θ
% where θ is angle between vectors, σ₁ and σ₂ are measurement uncertainties
% Simplified: σ_TRIAD ≈ σ_centroid / (sin(θ) * sqrt(2))
sigma_triad_factor = 1 / (sin(angle_sep_rad) * sqrt(2));
sigma_triad_arcsec = sigma_centroid_arcsec * sigma_triad_factor;

fprintf('  TRIAD CRLB:        %.2f arcsec (2 vectors, θ = %.1f°)\n', ...
    sigma_triad_arcsec, triad_info.angle_between_vectors);
fprintf('  Primary weight:    %.4f (most accurate vector)\n', triad_info.vector1_weight);
fprintf('  Secondary weight:  %.4f (ratio = %.2f:1)\n', ...
    triad_info.vector2_weight, weight_ratio);

% Check if actual error is close to theoretical limit
efficiency = sigma_triad_arcsec / angle_error_arcsec * 100;
fprintf('  Actual error:      %.2f arcsec\n', angle_error_arcsec);
fprintf('  Efficiency:        %.1f%% ', efficiency);
if efficiency > 80
    fprintf('✓ Near-optimal\n');
elseif efficiency > 50
    fprintf('○ Good\n');
else
    fprintf('⚠ Suboptimal (noise sources present)\n');
end

% TRIAD optimality assessment
fprintf('\n--- TRIAD Optimality Assessment ---\n');
if weight_ratio > 3.0
    fprintf('  Weight distribution: ✓ TRIAD optimal (σ₁ << σ₂, ratio %.1f:1)\n', weight_ratio);
    fprintf('  Primary vector dominates accuracy as expected.\n');
elseif weight_ratio > 1.5
    fprintf('  Weight distribution: ○ TRIAD acceptable (moderate ratio %.1f:1)\n', weight_ratio);
    fprintf('  QUEST may provide similar or slightly better accuracy.\n');
else
    fprintf('  Weight distribution: ⚠ Vectors have similar weights (%.1f:1)\n', weight_ratio);
    fprintf('  QUEST recommended: Better utilizes all %d observations.\n', triad_info.N_obs);
end

if triad_info.angle_between_vectors > 30 && triad_info.angle_between_vectors < 150
    fprintf('  Vector separation:   ✓ Good conditioning (%.1f°)\n', triad_info.angle_between_vectors);
elseif triad_info.angle_between_vectors > 15 && triad_info.angle_between_vectors < 165
    fprintf('  Vector separation:   ○ Acceptable (%.1f°)\n', triad_info.angle_between_vectors);
else
    fprintf('  Vector separation:   ⚠ Poor conditioning (%.1f°, near-collinear)\n', ...
        triad_info.angle_between_vectors);
end

if triad_info.N_obs > 2
    fprintf('  Additional stars:    %d observations not used by TRIAD\n', ...
        triad_info.N_obs - 2);
    fprintf('  Consider QUEST to utilize all available measurements.\n');
end

% Identify dominant noise source
if angle_error_arcsec < 2 * sigma_triad_arcsec
    dominant_noise = 'Centroid noise (measurement-limited)';
elseif triad_info.residual_mean > 100
    dominant_noise = 'Motion blur (dynamics-limited)';
elseif triad_info.residual_std / triad_info.residual_mean > 0.8
    dominant_noise = 'Outliers or inconsistent measurements';
else
    dominant_noise = 'Mixed noise sources';
end
fprintf('  Dominant source:     %s\n', dominant_noise);

% DCM Error Matrix
fprintf('\n--- Direction Cosine Matrix Error ---\n');
DCM_error = DCM_estimated * DCM_true'; % Should be close to identity
DCM_error_angles = [asin(DCM_error(3,2) - DCM_error(2,3)) / 2;
                    asin(DCM_error(1,3) - DCM_error(3,1)) / 2;
                    asin(DCM_error(2,1) - DCM_error(1,2)) / 2];
DCM_error_angles_arcsec = rad2deg(DCM_error_angles) * 3600;
fprintf('  Small-angle errors: [%.2f, %.2f, %.2f] arcsec\n', DCM_error_angles_arcsec);
fprintf('  Frobenius norm:     %.2e ', norm(DCM_error - eye(3), 'fro'));
if norm(DCM_error - eye(3), 'fro') < 1e-5
    fprintf('✓ Excellent\n');
elseif norm(DCM_error - eye(3), 'fro') < 1e-3
    fprintf('✓ Good\n');
else
    fprintf('○ Acceptable\n');
end

%% ========================================================================
%  7. RESULTS VISUALIZATION
% =========================================================================

% Generate comprehensive validation plots
saveFlag = 1;
plotTRIADResults(meas, N_STR, STR, DCM_true, DCM_estimated, ...
                 q_true, q_estimated, triad_info, saveFlag);

fprintf('\n=== All tasks completed successfully ===\n');