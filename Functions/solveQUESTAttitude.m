function [q_estimated, DCM_estimated, quest_info] = solveQUESTAttitude(meas, N_STR)
%==========================================================================
% solveQUESTAttitude: Solve Wahba's problem using QUEST algorithm for
%                     optimal attitude quaternion estimation.
%
% Inputs:
%   meas          - Structure array (1xN_STR) with star tracker measurements:
%                   .N_stars    - Number of detected stars
%                   .r_body     - 3xM unit vectors in body frame
%                   .r_ECI_ref  - 3xM reference unit vectors in ECI frame
%                   .weights    - Mx1 measurement weights (normalized)
%   N_STR         - Number of active star trackers.
%
% Outputs:
%   q_estimated   - Optimal attitude quaternion [qw; qx; qy; qz] (4x1).
%   DCM_estimated - Direction Cosine Matrix (ECI-to-Body, 3x3).
%   quest_info    - Structure with diagnostic information:
%                   .N_obs      - Total number of observations
%                   .lambda_max - Maximum eigenvalue of K-matrix
%                   .loss_wahba - Wahba loss function value
%                   .residuals  - Nx1 angular residuals [rad]
%                   .status     - String indicating solution quality
%
% Method:
%   1. Aggregate measurements from all STRs into global reference frames.
%   2. Construct attitude profile matrix B (weighted sum of dyadic products).
%   3. Build Davenport's K-matrix from B (4x4 symmetric).
%   4. Solve eigenvalue problem: K*q = λ*q for maximum eigenvalue.
%   5. Normalize and enforce positive scalar convention.
%
% References:
%   - Shuster, M. D. (1993). 
%     "A Survey of Attitude Representations". 
%     Navigation, 8(9), 439-517.
%   - Markley, F. L., & Mortari, D. (1999). 
%     "How to Estimate Attitude from Vector Observations". 
%     In Astrodynamics Specialist.
%
% Notes:
%   - Minimum 3 observations required (N_obs ≥ 3).
%   - Weights must sum to 1 (enforced internally).
%   - Scalar-first quaternion convention: [qw; qx; qy; qz].
%   - QUEST is algebraically optimal but sensitive to noise outliers.
%==========================================================================

    fprintf('\n=== Executing QUEST Algorithm ===\n');
    
    % Count total number of observations for pre-allocation
    N_total_obs = sum([meas.N_stars]);
    
    fprintf('  Total observations: %d stars\n', N_total_obs);
    
    % Validate sufficient observations
    if N_total_obs < 3
        error('Insufficient observations (%d < 3). QUEST requires ≥3 stars.', N_total_obs);
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
            
            r_body_all(:, idx_start:idx_end)  = meas(i_str).r_body;
            r_ECI_all(:, idx_start:idx_end)   = meas(i_str).r_ECI_ref;
            weights_all(idx_start:idx_end)    = meas(i_str).weights;
            
            idx_start = idx_end + 1;
        end
    end
    
    % Normalize weights to sum to unity
    weights_all = weights_all / sum(weights_all);
    
    % Construct attitude profile matrix B (Wahba's problem)
    % Convention: B = Σ w_i * (r_ECI * r_body'), for ECI-to-Body DCM
    B = zeros(3, 3);
    for i = 1:N_total_obs
        B = B + weights_all(i) * (r_ECI_all(:,i) * r_body_all(:,i)');
    end
    
    % Build Davenport's K-matrix (4x4 symmetric eigenvalue problem)
    S = B + B';                                  % Symmetric part
    Z = [B(2,3) - B(3,2);                        % Skew-symmetric extraction
         B(3,1) - B(1,3);
         B(1,2) - B(2,1)];
    sigma = trace(B);                            % Trace of B
    
    K = [S - sigma*eye(3),     Z; 
                       Z', sigma];
    
    % Solve eigenvalue problem for maximum eigenvalue
    [eigvec, eigval] = eig(K, 'vector');
    [lambda_max, idx_max] = max(eigval);
    
    % Extract optimal quaternion (eigenvector of λ_max)
    q_estimated = eigvec(:, idx_max);
    
    % Reorder to scalar-first convention: [qw; qx; qy; qz]
    q_estimated = [q_estimated(4); q_estimated(1:3)];
    
    % Normalize quaternion to unit magnitude
    q_estimated = q_estimated / norm(q_estimated);
    
    % Enforce positive scalar part convention
    if q_estimated(1) < 0
        q_estimated = -q_estimated;
    end
    
    fprintf('  Quaternion (est):   [%.6f, %.6f, %.6f, %.6f]\n', q_estimated);
    
    % Convert quaternion to DCM
    DCM_estimated = quat2dcm_custom(q_estimated);
    
    % Compute Wahba loss function and residuals
    loss_wahba = 0;
    residuals = zeros(N_total_obs, 1);
    for i = 1:N_total_obs
        residual_vec = r_body_all(:,i) - DCM_estimated * r_ECI_all(:,i);
        loss_wahba   = loss_wahba + weights_all(i) * norm(residual_vec)^2;
        residuals(i) = acos(max(min(dot(r_body_all(:,i), DCM_estimated * r_ECI_all(:,i)), 1), -1));
    end
    
    % Compute residual statistics in arcsec
    residual_mean_arcsec = rad2deg(mean(residuals)) * 3600;
    residual_std_arcsec  = rad2deg(std(residuals)) * 3600;
    residual_max_arcsec  = rad2deg(max(residuals)) * 3600;
    
    %% ====================================================================
    %  AUTOMATIC VALIDATION WITH STATUS INDICATORS
    % =====================================================================
    
    fprintf('\n--- Solution Quality Assessment ---\n');
    
    % 1. Lambda max validation (should be ≈ 1.0 after weight normalization)
    lambda_error = abs(lambda_max - 1.0);
    if lambda_error < 1e-6
        lambda_status = '✓ Excellent';
        lambda_color = '';
    elseif lambda_error < 1e-4
        lambda_status = '✓ Good';
        lambda_color = '';
    else
        lambda_status = '⚠ Suspicious';
        lambda_color = '(check weight normalization)';
    end
    fprintf('  Lambda (max):       %.10f   %s %s\n', lambda_max, lambda_status, lambda_color);
    
    % 2. Wahba loss validation (should be near-zero for good fit)
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
    fprintf('  Wahba loss:         %.8f     %s %s\n', loss_wahba, loss_status, loss_comment);
    
    % 3. Residual validation (depends on sensor specs and operational conditions)
    fprintf('  Residual (mean):     %.4f arcsec ', residual_mean_arcsec);
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
    
    fprintf('  Residual (std):     %.4f arcsec  ', residual_std_arcsec);
    if residual_std_arcsec < residual_mean_arcsec * 0.5
        fprintf('✓ Consistent\n');
    elseif residual_std_arcsec < residual_mean_arcsec * 1.0
        fprintf('○ Moderate scatter\n');
    else
        fprintf('⚠ High variance (possible outliers)\n');
    end
    
    fprintf('  Residual (max):     %.4f arcsec  ', residual_max_arcsec);
    if residual_max_arcsec < residual_mean_arcsec * 3
        fprintf('✓ No outliers\n');
    else
        fprintf('⚠ Potential outlier detected\n');
    end
    
    % 4. Quaternion norm validation
    q_norm = norm(q_estimated);
    fprintf('  Quaternion norm:    %.10f   ', q_norm);
    if abs(q_norm - 1.0) < 1e-10
        fprintf('✓ Unit quaternion\n');
    else
        fprintf('⚠ Non-unit (%.2e error)\n', abs(q_norm - 1.0));
    end
    
    % 5. DCM validation
    DCM_det = det(DCM_estimated);
    DCM_ortho_error = norm(DCM_estimated * DCM_estimated' - eye(3), 'fro');
    fprintf('  DCM determinant:    %.10f   ', DCM_det);
    if abs(DCM_det - 1.0) < 1e-10
        fprintf('✓ Valid rotation\n');
    else
        fprintf('⚠ Non-orthogonal (%.2e error)\n', abs(DCM_det - 1.0));
    end
    
    fprintf('  DCM orthogonality:  %.2e       ', DCM_ortho_error);
    if DCM_ortho_error < 1e-10
        fprintf('✓ Orthogonal\n');
    else
        fprintf('⚠ Error (check quat2dcm)\n');
    end
    
    %% ====================================================================
    %  PACKAGE DIAGNOSTIC INFORMATION
    % =====================================================================
    
    quest_info.N_obs            = N_total_obs;
    quest_info.lambda_max       = lambda_max;
    quest_info.loss_wahba       = loss_wahba;
    quest_info.residuals        = residuals;
    quest_info.residual_mean    = residual_mean_arcsec;
    quest_info.residual_std     = residual_std_arcsec;
    quest_info.residual_max     = residual_max_arcsec;
    quest_info.weights          = weights_all;
    quest_info.status           = overall_status;
    quest_info.q_norm           = q_norm;
    quest_info.DCM_determinant  = DCM_det;
    quest_info.DCM_ortho_error  = DCM_ortho_error;
    
end