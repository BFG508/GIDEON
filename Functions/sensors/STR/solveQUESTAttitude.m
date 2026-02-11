function [qEst, DCMEst, questInfo] = solveQUESTAttitude(meas, nSTR)
%==========================================================================
% solveQUESTAttitude: Solve Wahba's problem using QUEST algorithm for
%                     optimal attitude quaternion estimation.
%
% Inputs:
%   meas       - Structure array (1xN_STR) with star tracker measurements:
%                .nStars    - Number of detected stars
%                .rBody     - 3xM unit vectors in body frame
%                .rECI_ref  - 3xM reference unit vectors in ECI frame
%                .weights   - Mx1 measurement weights (normalized)
%   nSTR       - Number of active star trackers.
%
% Outputs:
%   qEst      - Optimal attitude quaternion [qw; qx; qy; qz] (4x1).
%   DCMEst    - Direction Cosine Matrix (ECI-to-Body, 3x3).
%   questInfo - Structure with diagnostic information:
%                .nObs      - Total number of observations
%                .lambdaMax - Maximum eigenvalue of K-matrix
%                .lossWahba - Wahba loss function value
%                .residuals - Nx1 angular residuals [rad]
%                .status    - String indicating solution quality
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
%   - Minimum 3 observations required (nObs ≥ 3).
%   - Weights must sum to 1 (enforced internally).
%   - Scalar-first quaternion convention: [qw; qx; qy; qz].
%   - QUEST is algebraically optimal but sensitive to noise outliers.
%==========================================================================

    fprintf('\n=== Executing QUEST Algorithm ===\n');
    
    % Count total number of observations for pre-allocation
    nTotalObs = sum([meas.nStars]);
    
    fprintf('  Total observations: %d stars\n', nTotalObs);
    
    % Validate sufficient observations
    if nTotalObs < 3
        error('Insufficient observations (%d < 3). QUEST requires ≥3 stars.', nTotalObs);
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
    
    % Construct attitude profile matrix B (Wahba's problem)
    % Convention: B = Σ w_i * (rECI * rBody'), for ECI-to-Body DCM
    B = zeros(3, 3);
    for i = 1:nTotalObs
        B = B + weights_all(i) * (rECI_all(:,i) * rBody_all(:,i)');
    end
    
    % Build Davenport's K-matrix (4x4 symmetric eigenvalue problem)
    S = B + B';                     % Symmetric part
    Z = [B(2,3) - B(3,2);           % Skew-symmetric extraction
         B(3,1) - B(1,3);
         B(1,2) - B(2,1)];
    sigma = trace(B);               % Trace of B
    
    K = [S - sigma*eye(3),     Z; 
                       Z', sigma];
    
    % Solve eigenvalue problem for maximum eigenvalue
    [eigVec, eigVal] = eig(K, 'vector');
    [lambdaMax, idxMax] = max(eigVal);
    
    % Extract optimal quaternion (eigenvector of λ_max)
    qEst = eigVec(:, idxMax);
    
    % Reorder to scalar-first convention: [qw; qx; qy; qz]
    qEst = [qEst(4); qEst(1:3)];
    
    % Normalize quaternion to unit magnitude
    qEst = qEst / norm(qEst);
    
    % Enforce positive scalar part convention
    if qEst(1) < 0
        qEst = - qEst;
    end
    
    fprintf('  Quaternion (est):   [%.6f, %.6f, %.6f, %.6f]\n', qEst);
    
    % Convert quaternion to DCM
    DCMEst = quat2dcm(qEst);
    
    % Compute Wahba loss function and residuals
    lossWahba = 0;
    residuals = zeros(nTotalObs, 1);
    for i = 1:nTotalObs
        residualVec  = rBody_all(:,i) - DCMEst * rECI_all(:,i);
        lossWahba    = lossWahba + weights_all(i) * norm(residualVec)^2;
        residuals(i) = acos(max(min(dot(rBody_all(:, i), DCMEst * rECI_all(:, i)), 1), -1));
    end
    
    % Compute residual statistics in arcsec
    residualMean = rad2deg(mean(residuals)) * 3600;
    residualStd  = rad2deg(std(residuals)) * 3600;
    residualMax  = rad2deg(max(residuals)) * 3600;
    
    %% ====================================================================
    %  AUTOMATIC VALIDATION WITH STATUS INDICATORS
    % =====================================================================
    
    fprintf('\n--- Solution Quality Assessment ---\n');
    
    % 1. Lambda max validation (should be ≈ 1.0 after weight normalization)
    lambdaErr = abs(lambdaMax - 1.0);
    if lambdaErr < 1e-6
        lambdaStatus = '✓ Excellent';
        lambdaColor = '';
    elseif lambdaErr < 1e-4
        lambdaStatus = '✓ Good';
        lambdaColor = '';
    else
        lambdaStatus = '⚠ Suspicious';
        lambdaColor = '(check weight normalization)';
    end
    fprintf('  Lambda (max):       %.10f   %s %s\n', lambdaMax, lambdaStatus, lambdaColor);
    
    % 2. Wahba loss validation (should be near-zero for good fit)
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
    fprintf('  Wahba loss:         %.8f     %s %s\n', lossWahba, lossStatus, lossComment);
    
    % 3. Residual validation (depends on sensor specs and operational conditions)
    fprintf('  Residual (mean):    %.4f arcsec ', residualMean);
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
    
    fprintf('  Residual (std):     %.4f arcsec  ', residualStd);
    if residualStd < residualMean * 0.5
        fprintf('✓ Consistent\n');
    elseif residualStd < residualMean * 1.0
        fprintf('○ Moderate scatter\n');
    else
        fprintf('⚠ High variance (possible outliers)\n');
    end
    
    fprintf('  Residual (max):     %.4f arcsec  ', residualMax);
    if residualMax < residualMean * 3
        fprintf('✓ No outliers\n');
    else
        fprintf('⚠ Potential outlier detected\n');
    end
    
    % 4. Quaternion norm validation
    qNorm = norm(qEst);
    fprintf('  Quaternion norm:    %.10f   ', qNorm);
    if abs(qNorm - 1.0) < 1e-10
        fprintf('✓ Unit quaternion\n');
    else
        fprintf('⚠ Non-unit (%.2e error)\n', abs(qNorm - 1.0));
    end
    
    % 5. DCM validation
    DCM_det = det(DCMEst);
    DCM_orthoError = norm(DCMEst * DCMEst' - eye(3), 'fro');
    fprintf('  DCM determinant:    %.10f   ', DCM_det);
    if abs(DCM_det - 1.0) < 1e-10
        fprintf('✓ Valid rotation\n');
    else
        fprintf('⚠ Non-orthogonal (%.2e error)\n', abs(DCM_det - 1.0));
    end
    
    fprintf('  DCM orthogonality:  %.2e       ', DCM_orthoError);
    if DCM_orthoError < 1e-10
        fprintf('✓ Orthogonal\n');
    else
        fprintf('⚠ Error (check quat2dcm)\n');
    end
    
    
    %% ====================================================================
    %  PACKAGE DIAGNOSTIC INFORMATION
    % =====================================================================
    
    questInfo.nObs            = nTotalObs;
    questInfo.lambdaMax       = lambdaMax;
    questInfo.lossWahba       = lossWahba;
    questInfo.residuals       = residuals;
    questInfo.residualMean    = residualMean;
    questInfo.residualStd     = residualStd;
    questInfo.residualMax     = residualMax;
    questInfo.weights         = weights_all;
    questInfo.status          = overallStatus;
    questInfo.qNorm           = qNorm;
    questInfo.DCM_determinant = DCM_det;
    questInfo.DCM_orthoErr    = DCM_orthoError;
end