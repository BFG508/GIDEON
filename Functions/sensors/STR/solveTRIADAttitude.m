function [qEst, DCMEst, triadInfo] = solveTRIADAttitude(meas, nSTR)
%==========================================================================
% solveTRIADAttitude: Solve attitude determination using TRIAD algorithm
%                     for two-vector observations.
%
% Inputs:
%   meas       - Structure array (1 x nSTR) with star tracker measurements:
%                .nStars        - Number of detected stars
%                .rBody         - 3xM unit vectors in body frame
%                .rECI_ref      - 3xM reference unit vectors in ECI frame
%                .weights        - Mx1 measurement weights (normalized)
%   nSTR       - Number of active star trackers.
%
% Outputs:
%   qEst      - Optimal attitude quaternion [qw; qx; qy; qz] (4x1).
%   DCMEst    - Direction Cosine Matrix (ECI-to-Body, 3x3).
%   triadInfo - Structure with diagnostic information:
%                .nObs          - Total number of observations
%                .lossWahba     - Wahba loss function value
%                .residuals      - Nx1 angular residuals [rad]
%                .status         - String indicating solution quality
%                .vector1_weight - Weight of primary observation
%                .vector2_weight - Weight of secondary observation
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
%   - Requires exactly 2 observations (nObs ≥ 2).
%   - Primary vector preserved exactly; secondary projected orthogonally.
%   - Scalar-first quaternion convention: [qw; qx; qy; qz].
%   - TRIAD is optimal when σ₁ << σ₂ (primary much more accurate).
%   - For equal-weighted observations, QUEST is generally superior.
%   - Observations must be non-collinear (v1 × v2 ≠ 0).
%==========================================================================

    fprintf('\n=== Executing TRIAD Algorithm ===\n');
    
    % Count total number of observations for validation
    nTotalObs = sum([meas.nStars]);
    fprintf(' Total observations: %d stars\n', nTotalObs);
    
    % Validate sufficient observations for TRIAD (minimum 2 required)
    if nTotalObs < 2
        error('Insufficient observations (%d < 2). TRIAD requires ≥ 2 stars.', nTotalObs);
    end
    
    % Pre-allocate arrays (avoids dynamic growth warning)
    rBody_all   = zeros(3, nTotalObs);
    rECI_all    = zeros(3, nTotalObs);
    weights_all = zeros(1, nTotalObs);
    
    % Aggregate measurements from all active star trackers
    idxStart = 1;
    for iSTR = 1:nSTR
        if meas(iSTR).nStars > 0
            idxEnd = idxStart + meas(iSTR).nStars - 1;
            
            rBody_all(:, idxStart:idxEnd) = meas(iSTR).rBody;
             rECI_all(:, idxStart:idxEnd) = meas(iSTR).rECI_ref;
             weights_all(idxStart:idxEnd) = meas(iSTR).weights;
            
            idxStart = idxEnd + 1;
        end
    end
    
    % Normalize weights to sum to unity
    weights_all = weights_all / sum(weights_all);
    
    % Select two highest-weighted observations (primary and secondary vectors)
    % Primary vector: most accurate measurement (highest weight)
    % Secondary vector: second most accurate measurement
    [weightsSorted, idxSorted] = sort(weights_all, 'descend');
    
    idxPrimary   = idxSorted(1); % Index of primary observation
    idxSecondary = idxSorted(2); % Index of secondary observation
    
    v1_ref  =  rECI_all(:, idxPrimary);   % Primary reference vector (ECI)
    v2_ref  =  rECI_all(:, idxSecondary); % Secondary reference vector (ECI)
    w1_body = rBody_all(:, idxPrimary);   % Primary body vector
    w2_body = rBody_all(:, idxSecondary); % Secondary body vector
    
    weightPrimary   = weightsSorted(1);
    weightSecondary = weightsSorted(2);
    
    fprintf(' Primary vector:   weight = %.4f (highest accuracy)\n', weightPrimary);
    fprintf(' Secondary vector: weight = %.4f\n'                   , weightSecondary);
    
    % Validate non-collinearity (vectors must not be parallel)
    crossRef       = cross(v1_ref, v2_ref);
    crossBody      = cross(w1_body, w2_body);
    crossRef_norm  = norm(crossRef);
    crossBody_norm = norm(crossBody);
    
    if crossRef_norm < 1e-6 || crossBody_norm < 1e-6
        error(['Collinear vectors detected (||v1 × v2|| = %.2e). ', ...
               'TRIAD requires non-parallel observations.'], ...
               min(crossRef_norm, crossBody_norm));
    end
    
    % Compute angle between vectors (for diagnostics)
    angleBetween = acosd(max(min(dot(v1_ref, v2_ref), 1), -1));
    fprintf(' Angle between vectors: %.2f deg ', angleBetween);
    if angleBetween > 30 && angleBetween < 150
        fprintf('✓ Good separation\n');
    elseif angleBetween > 15 && angleBetween < 165
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
    
    r1_ref = v1_ref / norm(v1_ref);    % Normalize primary (redundant if already unit)
    r2_ref = crossRef / crossRef_norm; % Perpendicular direction
    r3_ref = cross(r1_ref, r2_ref);    % Complete triad
    
    % Body frame triad (satellite frame):
    %   w1 = w1_body / ||w1_body||    (primary direction, preserved exactly)
    %   w2 = (w1 × w2) / ||w1 × w2||  (perpendicular to plane of w1 and w2)
    %   w3 = w1 × w2                  (completes right-handed orthonormal triad)
    
    w1_bodyNorm = w1_body / norm(w1_body);         % Normalize primary (redundant if already unit)
    w2_bodyNorm = crossBody / crossBody_norm;      % Perpendicular direction
    w3_bodyNorm = cross(w1_bodyNorm, w2_bodyNorm); % Complete triad
    
    % Assemble triads as column matrices
    M_ref  = [r1_ref, r2_ref, r3_ref];                % Reference triad (3×3)
    M_body = [w1_bodyNorm, w2_bodyNorm, w3_bodyNorm]; % Body triad (3×3)
    
    %% ====================================================================
    % COMPUTE ATTITUDE MATRIX (Direction Cosine Matrix)
    % =====================================================================
    % DCM transforms reference frame to body frame:
    %   w_i = A * r_i  for i = 1, 2, 3
    % Therefore:
    %   A = [w1 w2 w3] * [r1 r2 r3]ᵀ = M_body * M_refᵀ
    
    DCMEst = M_body * M_ref';
    
    %% ====================================================================
    % EXTRACT QUATERNION FROM DCM (Shepperd's Method)
    % =====================================================================
    % Use custom DCM-to-quaternion conversion (already available in codebase)
    qEst = dcm2quat(DCMEst);
    
    fprintf(' Quaternion (est): [%.6f, %.6f, %.6f, %.6f]\n', qEst);
    
    %% ====================================================================
    % COMPUTE WAHBA LOSS FUNCTION AND RESIDUALS
    % =====================================================================
    % Evaluate fit quality by computing residuals for all observations
    % (not just the two used in TRIAD construction)
    
    lossWahba = 0;
    residuals = zeros(nTotalObs, 1);
    
    for i = 1:nTotalObs
        % Residual vector: difference between observed and predicted body vector
        residual_vec = rBody_all(:,i) - DCMEst * rECI_all(:,i);
        lossWahba = lossWahba + weights_all(i) * norm(residual_vec)^2;
        
        % Angular residual: angle between observed and predicted directions
        cos_angle = dot(rBody_all(:,i), DCMEst * rECI_all(:,i));
        cos_angle = max(min(cos_angle, 1), -1); % Clamp to [-1, 1] for numerical stability
        residuals(i) = acos(cos_angle);
    end
    
    % Compute residual statistics in arcsec
    residualMean = rad2deg(mean(residuals)) * 3600;
    residualStd  = rad2deg(std(residuals)) * 3600;
    residualMax  = rad2deg(max(residuals)) * 3600;
    
    %% ====================================================================
    % AUTOMATIC VALIDATION WITH STATUS INDICATORS
    % =====================================================================
    fprintf('\n--- Solution Quality Assessment ---\n');
    
    % 1. Wahba loss validation (should be near-zero for good fit)
    if lossWahba < 1e-6
        lossStatus = '✓ Excellent';
        lossComment = '(near-perfect fit)';
    elseif lossWahba < 1e-3
        lossStatus = '✓ Good';
        lossComment = '(low noise)';
    elseif lossWahba < 0.1
        lossStatus = '○ Acceptable';
        lossComment = '(moderate noise/motion blur)';
    elseif lossWahba < 1.0
        lossStatus = '⚠ High';
        lossComment = '(significant noise/misalignment)';
    else
        lossStatus = '❌ Poor';
        lossComment = '(check measurement consistency)';
    end
    fprintf(' Wahba loss:        %.8f    %s %s\n', lossWahba, lossStatus, lossComment);
    
    % 2. Residual validation (depends on sensor specs and operational conditions)
    fprintf(' Residual (mean):   %.4f arcsec ', residualMean);
    if residualMean < 10
        fprintf('✓ Excellent (better than typical STR)\n');
        overallStatus = 'EXCELLENT';
    elseif residualMean < 30
        fprintf('✓ Good (within high-accuracy STR specs)\n');
        overallStatus = 'GOOD';
    elseif residualMean < 100
        fprintf('○ Acceptable (within standard STR specs)\n');
        overallStatus = 'ACCEPTABLE';
    elseif residualMean < 500
        fprintf('⚠ High (motion blur or moderate noise)\n');
        overallStatus = 'HIGH NOISE';
    else
        fprintf('❌ Poor (check data quality)\n');
        overallStatus = 'POOR';
    end
    
    fprintf(' Residual (std):    %.4f arcsec ', residualStd);
    if residualStd < residualMean * 0.5
        fprintf('✓ Consistent\n');
    elseif residualStd < residualMean * 1.0
        fprintf('○ Moderate scatter\n');
    else
        fprintf('⚠ High variance (possible outliers)\n');
    end
    
    fprintf(' Residual (max):    %.4f arcsec ', residualMax);
    if residualMax < residualMean * 3
        fprintf('✓ No outliers\n');
    else
        fprintf('⚠ Potential outlier detected\n');
    end
    
    % 3. Quaternion norm validation
    qNorm = norm(qEst);
    fprintf(' Quaternion norm:   %.10f  ', qNorm);
    if abs(qNorm - 1.0) < 1e-10
        fprintf('✓ Unit quaternion\n');
    else
        fprintf('⚠ Non-unit (%.2e error)\n', abs(qNorm - 1.0));
    end
    
    % 4. DCM validation (orthogonality and determinant)
    DCM_det = det(DCMEst);
    DCM_orthoErr = norm(DCMEst * DCMEst' - eye(3), 'fro');
    
    fprintf(' DCM determinant:   %.10f  ', DCM_det);
    if abs(DCM_det - 1.0) < 1e-10
        fprintf('✓ Valid rotation\n');
    else
        fprintf('⚠ Non-orthogonal (%.2e error)\n', abs(DCM_det - 1.0));
    end
    
    fprintf(' DCM orthogonality: %.2e      ', DCM_orthoErr);
    if DCM_orthoErr < 1e-10
        fprintf('✓ Orthogonal\n');
    else
        fprintf('⚠ Error (check triad construction)\n');
    end
    
    % 5. TRIAD-specific validation: Check if unused observations are consistent
    if nTotalObs > 2
        fprintf('\n--- Additional Observations Assessment ---\n');
        residualsUnused    = residuals(idxSorted(3:end));
        residualUnusedMean = rad2deg(mean(residualsUnused)) * 3600;
        fprintf(' Unused obs residual (mean): %.4f arcsec ', residualUnusedMean);
        if residualUnusedMean < 1.5 * residualMean
            fprintf('✓ Consistent with TRIAD solution\n');
        else
            fprintf('⚠ Higher than primary/secondary (consider QUEST)\n');
        end
    end
    
    %% ====================================================================
    % PACKAGE DIAGNOSTIC INFORMATION
    % =====================================================================

    triadInfo.nObs            = nTotalObs;
    triadInfo.lossWahba       = lossWahba;
    triadInfo.residuals       = residuals;
    triadInfo.residualMean    = residualMean;
    triadInfo.residualStd     = residualStd;
    triadInfo.residualMax     = residualMax;
    triadInfo.weights         = weights_all;
    triadInfo.status          = overallStatus;
    triadInfo.qNorm           = qNorm;
    triadInfo.DCM_determinant = DCM_det;
    triadInfo.DCM_orthoErr    = DCM_orthoErr;
    triadInfo.weightPrimary   = weightPrimary;
    triadInfo.weightSecondary = weightSecondary;
    triadInfo.angleBetween    = angleBetween;
    triadInfo.idxPrimary      = idxPrimary;
    triadInfo.idxSecondary    = idxSecondary;
end