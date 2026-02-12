function EKF = updateSTR(EKF, strMeas, nSTR)
%==========================================================================
% updateSTR: Star Tracker measurement update step for MEKF. Processes matched
%            star pairs from all STR units, builds stacked measurement model,
%            and performs batch Kalman update using cross-product residuals.
%
% INPUTS:
%    EKF     - EKF structure with fields:
%        .qNom  - Nominal attitude quaternion (4x1), ECI to Body
%        .x     - Error state vector (6x1): [deltaTheta; biasGyro]
%        .P     - Error state covariance matrix (6x6)
%        .R_str - Star Tracker measurement noise covariance (3x3) [rad²]
%
%    strMeas - Cell array of STR measurements (1 x nSTR)
%              Each cell contains struct with fields:
%        .nStars - Array of matched star pairs with:
%            .catVec - Catalog direction in ECI frame (3x1, unit vector)
%            .bVec   - Measured direction in body frame (3x1, unit vector)
%
%    nSTR    - Number of star tracker units
%
% OUTPUTS:
%    EKF     - Updated EKF structure with corrected:
%        .qNom - Nominal quaternion after attitude correction
%        .x    - Error state (attitude error reset to zero)
%        .P    - Updated error covariance
%
% METHOD:
%    1. Collect all valid star measurements from all STR units
%    2. For each star: compute cross-product residual y = s_meas × s_pred
%    3. Project 3D residuals to 2D tangent plane (2 measurements per star)
%    4. Build stacked measurement model: y_stack = H_stack * x + noise
%    5. Batch Kalman update with all stars simultaneously
%    6. Apply correction to quaternion and reset attitude error
%
% REFERENCE:
%    Markley & Crassidis, "Fundamentals of Spacecraft Attitude 
%    Determination and Control", Section 7.6
%
%==========================================================================

    % --- 1. Collect All Valid Star Measurements ---
    allStarsCell = cell(nSTR, 1);
    totalStars   = 0;
    
    for iSTR = 1:nSTR
        if ~isempty(strMeas(iSTR)) && isfield(strMeas(iSTR), 'nStars')
            allStarsCell{iSTR} = strMeas(iSTR).nStars;
            totalStars = totalStars + length(strMeas(iSTR).nStars);
        end
    end
    
    if totalStars < 3
        return;
    end
    
    % Concatenate all stars
    allStars = vertcat(allStarsCell{:});
    nStars   = length(allStars);
    
    % --- 2. Build Stacked Measurement Vectors ---
    nMeas   = nStars * 2;
    yStack  = zeros(nMeas, 1);
    HStack  = zeros(nMeas, 6);
    DCM_nom = quat2dcm(EKF.qNom);
    
    for i = 1:nStars
        % Catalog and measured star directions
        sCat  = allStars(i).catVec;
        sMeas = allStars(i).bVec;
        
        % Predicted measurement
        sPred = DCM_nom * sCat;
        
        % Innovation (cross-product residual)
        yFull = cross(sMeas, sPred);
        
        % Project to 2D tangent plane
        [y2D, T] = projectToTangPlane(yFull, sPred);
        
        % Measurement Jacobian
        H3D = computeMeasJacobian_Att(sPred);
        H2D = [T * H3D(:,1:3), zeros(2,3)];
        
        % Stack measurements
        idx = (i-1)*2 + (1:2);
        yStack(idx)   = y2D;
        HStack(idx,:) = H2D;
    end
    
    % --- 3. Build Stacked Measurement Noise Covariance ---
    RStack = kron(eye(nStars), EKF.R_str(1:2,1:2));
    
    % --- 4. Kalman Gain and State Update ---
    S     = HStack * EKF.P * HStack' + RStack;
    K     = EKF.P * HStack' / S;
    EKF.x = EKF.x + K * yStack;
    
    % --- 5. Update Covariance (Joseph Form) ---
    I     = eye(6);
    IKH   = I - K * HStack;
    EKF.P = IKH * EKF.P * IKH' + K * RStack * K';
    
    % --- 6. Apply Correction to Nominal Quaternion ---
    deltaQ   = smallAng2quat(EKF.x(1:3));
    EKF.qNom = quatmultiply(deltaQ, EKF.qNom);
    EKF.qNom = EKF.qNom / norm(EKF.qNom);
    
    % --- 7. Reset Attitude Error (MEKF Property) ---
    EKF.x(1:3) = zeros(3,1);
end

function [y2D, T] = projectToTangPlane(yFull, sPred)
%==========================================================================
% projectToTangPlane: Projects 3D cross-product residual onto 2D tangent
%                     plane perpendicular to predicted star direction.
%                     Extracts two independent measurements per star.
%
% INPUTS:
%    yFull - 3D cross-product residual (3x1)
%    sPred - Predicted star direction (3x1, unit vector)
%            Defines normal to tangent plane
%
% OUTPUTS:
%    y2D   - 2D projected residual (2x1)
%    T     - Projection matrix (2x3)
%==========================================================================

    % --- 1. Choose First Tangent Vector ---
    if abs(sPred(1)) < 0.9
        v = [1; 0; 0];
    else
        v = [0; 1; 0];
    end
    
    % --- 2. Gram-Schmidt Orthogonalization ---
    e1 = v - (v' * sPred) * sPred;
    e1 = e1 / norm(e1);
    
    % --- 3. Second Tangent Vector (Cross Product) ---
    e2 = cross(sPred, e1);
    e2 = e2 / norm(e2);
    
    % --- 4. Projection Matrix and Residual ---
    T   = [e1'; e2'];
    y2D = T * yFull;

end