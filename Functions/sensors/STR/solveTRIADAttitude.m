function [q_estimated, DCM_estimated, triad_info] = solveTRIADAttitude(meas, N_STR)
%==========================================================================
% solveTRIADAttitude: Solve attitude determination using TRIAD algorithm
%                     for two-vector observations.
%
% Inputs:
%   meas          - Structure array (1xN_STR) with star tracker measurements:
%                   .N_stars        - Number of detected stars
%                   .r_body         - 3xM unit vectors in body frame
%                   .r_ECI_ref      - 3xM reference unit vectors in ECI frame
%                   .weights        - Mx1 measurement weights (normalized)
%   N_STR         - Number of active star trackers.
%
% Outputs:
%   q_estimated   - Optimal attitude quaternion [qw; qx; qy; qz] (4x1).
%   DCM_estimated - Direction Cosine Matrix (ECI-to-Body, 3x3).
%   triad_info    - Structure with diagnostic information:
%                   .N_obs          - Total number of observations
%                   .loss_wahba     - Wahba loss function value
%                   .residuals      - Nx1 angular residuals [rad]
%                   .status         - String indicating solution quality
%                   .vector1_weight - Weight of primary observation
%                   .vector2_weight - Weight of secondary observation
%
% Method:
%   1. Select two highest-weighted observations (primary and secondary).
%   2. Construct orthonormal triads in both reference and body frames:
%      - r1 = v1 (primary reference vector, normalized)
%      - r2 = (v1 × v2) / ||v1 × v2|| (perpendicular direction)
%      - r3 = r1 × r2 (completes right-handed triad)
%      - Same construction for body frame: w1, w2, w3
%   3. Compute DCM directly: A = [w1 w2 w3] * [r1 r2 r3]ᵀ
%   4. Extract quaternion from DCM using Shepperd's method.
%
% References:
%   - Black, H. D. (1964). 
%     "A Passive System for Determining the Attitude of a Satellite". 
%     AIAA journal, 2(7), 1350-1351.
%   - Shuster, M. D., & Oh, S. D. (1981). 
%     "Three-Axis Attitude Determination from Vector Observations". 
%     Journal of guidance and Control, 4(1), 70-77.
%   - Shuster, M. D. (2007). 
%     "The Optimization of TRIAD". 
%     The Journal of the Astronautical Sciences, 55(2), 245-257.
%
% Notes:
%   - Requires exactly 2 observations (N_obs ≥ 2).
%   - Primary vector preserved exactly; secondary projected orthogonally.
%   - Scalar-first quaternion convention: [qw; qx; qy; qz].
%   - TRIAD is optimal when σ₁ << σ₂ (primary much more accurate).
%   - For equal-weighted observations, QUEST is generally superior.
%   - Observations must be non-collinear (v1 × v2 ≠ 0).
%==========================================================================

    fprintf('\n=== Executing TRIAD Algorithm ===\n');
    
    % Count total number of observations for validation
    N_total_obs = sum([meas.N_stars]);
    fprintf(' Total observations: %d stars\n', N_total_obs);
    
    % Validate sufficient observations for TRIAD (minimum 2 required)
    if N_total_obs < 2
        error('Insufficient observations (%d < 2). TRIAD requires ≥ 2 stars.', N_total_obs);
    end
    
    % Pre-allocate arrays (avoids dynamic growth warning)
    r_body_all  = zeros(3, N_total_obs);
    r_ECI_all   = zeros(3, N_total_obs);
    weights_all = zeros(1, N_total_obs);
    
    % Aggregate measurements from all active star trackers
    idx_start = 1;
    for i_str = 1:N_STR
        if meas(i_str).N_stars > 0
            idx_end = idx_start + meas(i_str).N_stars - 1;
            r_body_all(:, idx_start:idx_end) = meas(i_str).r_body;
            r_ECI_all(:, idx_start:idx_end) = meas(i_str).r_ECI_ref;
            weights_all(idx_start:idx_end) = meas(i_str).weights;
            idx_start = idx_end + 1;
        end
    end
    
    % Normalize weights to sum to unity
    weights_all = weights_all / sum(weights_all);
    
    % Select two highest-weighted observations (primary and secondary vectors)
    % Primary vector: most accurate measurement (highest weight)
    % Secondary vector: second most accurate measurement
    [weights_sorted, idx_sorted] = sort(weights_all, 'descend');
    
    idx_primary   = idx_sorted(1); % Index of primary observation
    idx_secondary = idx_sorted(2); % Index of secondary observation
    
    v1_ref  =  r_ECI_all(:, idx_primary);   % Primary reference vector (ECI)
    v2_ref  =  r_ECI_all(:, idx_secondary); % Secondary reference vector (ECI)
    w1_body = r_body_all(:, idx_primary);   % Primary body vector
    w2_body = r_body_all(:, idx_secondary); % Secondary body vector
    
    weight_primary = weights_sorted(1);
    weight_secondary = weights_sorted(2);
    
    fprintf(' Primary vector:   weight = %.4f (highest accuracy)\n', weight_primary);
    fprintf(' Secondary vector: weight = %.4f\n', weight_secondary);
    
    % Validate non-collinearity (vectors must not be parallel)
    cross_ref       = cross(v1_ref, v2_ref);
    cross_body      = cross(w1_body, w2_body);
    cross_ref_norm  = norm(cross_ref);
    cross_body_norm = norm(cross_body);
    
    if cross_ref_norm < 1e-6 || cross_body_norm < 1e-6
        error(['Collinear vectors detected (||v1 × v2|| = %.2e). ', ...
               'TRIAD requires non-parallel observations.'], ...
               min(cross_ref_norm, cross_body_norm));
    end
    
    % Compute angle between vectors (for diagnostics)
    angle_between = acos(max(min(dot(v1_ref, v2_ref), 1), -1));
    angle_between_deg = rad2deg(angle_between);
    fprintf(' Angle between vectors: %.2f deg ', angle_between_deg);
    if angle_between_deg > 30 && angle_between_deg < 150
        fprintf('✓ Good separation\n');
    elseif angle_between_deg > 15 && angle_between_deg < 165
        fprintf('○ Acceptable separation\n');
    else
        fprintf('⚠ Poor separation (near-collinear)\n');
    end
    
    %% ====================================================================
    % CONSTRUCT ORTHONORMAL TRIADS (TRIAD Algorithm Core)
    % =====================================================================
    % Reference frame triad (ECI frame):
    %   r1 = v1 / ||v1||              (primary direction, preserved exactly)
    %   r2 = (v1 × v2) / ||v1 × v2||  (perpendicular to plane of v1 and v2)
    %   r3 = r1 × r2                  (completes right-handed orthonormal triad)
    
    r1_ref = v1_ref / norm(v1_ref);      % Normalize primary (redundant if already unit)
    r2_ref = cross_ref / cross_ref_norm; % Perpendicular direction
    r3_ref = cross(r1_ref, r2_ref);      % Complete triad
    
    % Body frame triad (satellite frame):
    %   w1 = w1_body / ||w1_body||    (primary direction, preserved exactly)
    %   w2 = (w1 × w2) / ||w1 × w2||  (perpendicular to plane of w1 and w2)
    %   w3 = w1 × w2                  (completes right-handed orthonormal triad)
    
    w1_body_norm = w1_body / norm(w1_body);           % Normalize primary (redundant if already unit)
    w2_body_norm = cross_body / cross_body_norm;      % Perpendicular direction
    w3_body_norm = cross(w1_body_norm, w2_body_norm); % Complete triad
    
    % Assemble triads as column matrices
    M_ref  = [r1_ref, r2_ref, r3_ref];                   % Reference triad (3×3)
    M_body = [w1_body_norm, w2_body_norm, w3_body_norm]; % Body triad (3×3)
    
    %% ====================================================================
    % COMPUTE ATTITUDE MATRIX (Direction Cosine Matrix)
    % =====================================================================
    % DCM transforms reference frame to body frame:
    %   w_i = A * r_i  for i = 1, 2, 3
    % Therefore:
    %   A = [w1 w2 w3] * [r1 r2 r3]ᵀ = M_body * M_refᵀ
    
    DCM_estimated = M_body * M_ref';
    
    %% ====================================================================
    % EXTRACT QUATERNION FROM DCM (Shepperd's Method)
    % =====================================================================
    % Use custom DCM-to-quaternion conversion (already available in codebase)
    q_estimated = dcm2quat(DCM_estimated);
    
    fprintf(' Quaternion (est): [%.6f, %.6f, %.6f, %.6f]\n', q_estimated);
    
    %% ====================================================================
    % COMPUTE WAHBA LOSS FUNCTION AND RESIDUALS
    % =====================================================================
    % Evaluate fit quality by computing residuals for all observations
    % (not just the two used in TRIAD construction)
    
    loss_wahba = 0;
    residuals = zeros(N_total_obs, 1);
    
    for i = 1:N_total_obs
        % Residual vector: difference between observed and predicted body vector
        residual_vec = r_body_all(:,i) - DCM_estimated * r_ECI_all(:,i);
        loss_wahba = loss_wahba + weights_all(i) * norm(residual_vec)^2;
        
        % Angular residual: angle between observed and predicted directions
        cos_angle = dot(r_body_all(:,i), DCM_estimated * r_ECI_all(:,i));
        cos_angle = max(min(cos_angle, 1), -1); % Clamp to [-1, 1] for numerical stability
        residuals(i) = acos(cos_angle);
    end
    
    % Compute residual statistics in arcsec
    residual_mean_arcsec = rad2deg(mean(residuals)) * 3600;
    residual_std_arcsec = rad2deg(std(residuals)) * 3600;
    residual_max_arcsec = rad2deg(max(residuals)) * 3600;
    
    %% ====================================================================
    % AUTOMATIC VALIDATION WITH STATUS INDICATORS
    % =====================================================================
    fprintf('\n--- Solution Quality Assessment ---\n');
    
    % 1. Wahba loss validation (should be near-zero for good fit)
    if loss_wahba < 1e-6
        loss_status = '✓ Excellent';
        loss_comment = '(near-perfect fit)';
    elseif loss_wahba < 1e-3
        loss_status = '✓ Good';
        loss_comment = '(low noise)';
    elseif loss_wahba < 0.1
        loss_status = '○ Acceptable';
        loss_comment = '(moderate noise/motion blur)';
    elseif loss_wahba < 1.0
        loss_status = '⚠ High';
        loss_comment = '(significant noise/misalignment)';
    else
        loss_status = '❌ Poor';
        loss_comment = '(check measurement consistency)';
    end
    fprintf(' Wahba loss:      %.8f    %s %s\n', loss_wahba, loss_status, loss_comment);
    
    % 2. Residual validation (depends on sensor specs and operational conditions)
    fprintf(' Residual (mean): %.4f arcsec ', residual_mean_arcsec);
    if residual_mean_arcsec < 10
        fprintf('✓ Excellent (better than typical STR)\n');
        overall_status = 'EXCELLENT';
    elseif residual_mean_arcsec < 30
        fprintf('✓ Good (within high-accuracy STR specs)\n');
        overall_status = 'GOOD';
    elseif residual_mean_arcsec < 100
        fprintf('○ Acceptable (within standard STR specs)\n');
        overall_status = 'ACCEPTABLE';
    elseif residual_mean_arcsec < 500
        fprintf('⚠ High (motion blur or moderate noise)\n');
        overall_status = 'HIGH_NOISE';
    else
        fprintf('❌ Poor (check data quality)\n');
        overall_status = 'POOR';
    end
    
    fprintf(' Residual (std):  %.4f arcsec ', residual_std_arcsec);
    if residual_std_arcsec < residual_mean_arcsec * 0.5
        fprintf('✓ Consistent\n');
    elseif residual_std_arcsec < residual_mean_arcsec * 1.0
        fprintf('○ Moderate scatter\n');
    else
        fprintf('⚠ High variance (possible outliers)\n');
    end
    
    fprintf(' Residual (max):  %.4f arcsec ', residual_max_arcsec);
    if residual_max_arcsec < residual_mean_arcsec * 3
        fprintf('✓ No outliers\n');
    else
        fprintf('⚠ Potential outlier detected\n');
    end
    
    % 3. Quaternion norm validation
    q_norm = norm(q_estimated);
    fprintf(' Quaternion norm: %.10f  ', q_norm);
    if abs(q_norm - 1.0) < 1e-10
        fprintf('✓ Unit quaternion\n');
    else
        fprintf('⚠ Non-unit (%.2e error)\n', abs(q_norm - 1.0));
    end
    
    % 4. DCM validation (orthogonality and determinant)
    DCM_det = det(DCM_estimated);
    DCM_ortho_error = norm(DCM_estimated * DCM_estimated' - eye(3), 'fro');
    
    fprintf(' DCM determinant: %.10f  ', DCM_det);
    if abs(DCM_det - 1.0) < 1e-10
        fprintf('✓ Valid rotation\n');
    else
        fprintf('⚠ Non-orthogonal (%.2e error)\n', abs(DCM_det - 1.0));
    end
    
    fprintf(' DCM orthogonality: %.2e    ', DCM_ortho_error);
    if DCM_ortho_error < 1e-10
        fprintf('✓ Orthogonal\n');
    else
        fprintf('⚠ Error (check triad construction)\n');
    end
    
    % 5. TRIAD-specific validation: Check if unused observations are consistent
    if N_total_obs > 2
        fprintf('\n--- Additional Observations Assessment ---\n');
        residuals_unused = residuals(idx_sorted(3:end));
        residual_unused_mean = rad2deg(mean(residuals_unused)) * 3600;
        fprintf(' Unused obs residual (mean): %.4f arcsec ', residual_unused_mean);
        if residual_unused_mean < 1.5 * residual_mean_arcsec
            fprintf('✓ Consistent with TRIAD solution\n');
        else
            fprintf('⚠ Higher than primary/secondary (consider QUEST)\n');
        end
    end
    
    %% ====================================================================
    % PACKAGE DIAGNOSTIC INFORMATION
    % =====================================================================
    triad_info.N_obs = N_total_obs;
    triad_info.loss_wahba = loss_wahba;
    triad_info.residuals = residuals;
    triad_info.residual_mean = residual_mean_arcsec;
    triad_info.residual_std = residual_std_arcsec;
    triad_info.residual_max = residual_max_arcsec;
    triad_info.weights = weights_all;
    triad_info.status = overall_status;
    triad_info.q_norm = q_norm;
    triad_info.DCM_determinant = DCM_det;
    triad_info.DCM_ortho_error = DCM_ortho_error;
    triad_info.vector1_weight = weight_primary;
    triad_info.vector2_weight = weight_secondary;
    triad_info.angle_between_vectors = angle_between_deg;
    triad_info.idx_primary = idx_primary;
    triad_info.idx_secondary = idx_secondary;
end