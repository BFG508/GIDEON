function ekf = updateSTR(ekf, strMeas, nSTR)
%==========================================================================
% updateSTR - Star Tracker measurement update for MEKF
%
% INPUTS:
%   ekf     - EKF state structure (from initializeMEKF or predictMEKF)
%   strMeas - Cell array of STR measurements from all units
%             Each cell contains struct with fields:
%               .stars - Array of matched star pairs (catalog vs measured)
%   nSTR    - Number of star tracker units
%
% OUTPUTS:
%   ekf - Updated EKF state structure with corrected attitude and covariance
%
% ALGORITHM:
%   For each matched star pair (catalog direction vs. measured direction):
%     1. Predict measured direction using nominal attitude
%     2. Compute innovation (cross-product residual)
%     3. Build stacked measurement Jacobian
%     4. Perform batch Kalman update with all star measurements
%     5. Apply correction to nominal quaternion and reset error state
%
% NOTE: This function uses star vector pairs directly, which is more
%       robust than using STR-computed quaternions.
%==========================================================================

    %% 1. COLLECT ALL VALID STAR MEASUREMENTS
    allStars = [];
    
    for iSTR = 1:nSTR
        if ~isempty(strMeas(iSTR)) && isfield(strMeas(iSTR), 'nStars')
            % Concatenate stars from this STR unit
            allStars = [allStars; strMeas(iSTR).nStars];
        end
    end
    
    % Check if we have any valid measurements
    if isempty(allStars)
        % No STR data available, skip update
        return;
    end
    
    nStars = length(allStars);
    
    % Require minimum number of stars for robust update
    if nStars < 3
        % Not enough stars, skip update (need ≥3 for 3-axis attitude)
        return;
    end
    
    %% 2. BUILD STACKED MEASUREMENT VECTORS
    % Each star provides a 2D constraint (cross-product residual)
    % We use 2 axes of the 3D residual (the 3rd is zero by construction)
    
    nMeas = nStars * 2;  % 2 measurements per star
    y_stack = zeros(nMeas, 1);  % Stacked innovation
    H_stack = zeros(nMeas, 6);  % Stacked Jacobian
    
    DCM_nom = quat2dcm(ekf.q_nom);
    
    for i = 1:nStars
        % Catalog direction (ECI frame, unit vector)
        s_cat = allStars(i).catVec;  % 3x1
        
        % Measured direction (body frame, unit vector)
        s_meas = allStars(i).bVec;   % 3x1
        
        %% Predicted measurement
        % Rotate catalog direction to body frame using nominal attitude
        s_pred = DCM_nom * s_cat;  % 3x1
        
        %% Innovation (cross-product residual)
        % Full 3D residual: y = s_meas × s_pred
        % This captures the attitude error in the plane perpendicular to s_pred
        y_full = cross(s_meas, s_pred);  % 3x1
        
        % Extract 2D residual (discard component along s_pred, which is ~0)
        % Use a projection to get 2 independent measurements
        [y_2D, T] = projectToTangentPlane(y_full, s_pred);
        
        %% Measurement Jacobian
        % H = ∂y/∂δθ where y = s_meas × C(q_nom) * s_cat
        % Full 3D Jacobian: H_3D = -[s_pred]×
        % Projected to 2D: H_2D = T * H_3D
        H_3D = computeSTRJacobian(s_pred);  % 3x6
        H_2D = T * H_3D(:,1:3);             % 2x3 (only attitude part)
        H_2D = [H_2D, zeros(2,3)];          % 2x6 (add zero bias part)
        
        %% Stack measurements
        idx = (i-1)*2 + (1:2);
        y_stack(idx) = y_2D;
        H_stack(idx,:) = H_2D;
    end
    
    %% 3. BUILD STACKED MEASUREMENT NOISE COVARIANCE
    % Each star has independent noise, R_star is 2x2 per star
    % Stack into block-diagonal matrix
    R_stack = kron(eye(nStars), ekf.R_str(1:2,1:2));  % nMeas x nMeas
    
    %% 4. KALMAN GAIN AND STATE UPDATE
    S = H_stack * ekf.P * H_stack' + R_stack;  % Innovation covariance
    K = ekf.P * H_stack' / S;                  % Kalman gain, 6 x nMeas
    
    % Update error state
    ekf.x = ekf.x + K * y_stack;
    
    %% 5. UPDATE COVARIANCE (Joseph Form)
    I   = eye(6);
    IKH = I - K * H_stack;
    ekf.P = IKH * ekf.P * IKH' + K * R_stack * K';
    
    %% 6. RESET: Apply Correction to Nominal Quaternion
    delta_q = smallAngle2Quat(ekf.x(1:3));
    ekf.q_nom = quatmultiply(delta_q, ekf.q_nom);
    ekf.q_nom = ekf.q_nom / norm(ekf.q_nom);  % Renormalize
    
    % Reset attitude error
    ekf.x(1:3) = zeros(3,1);
end

%% =======================================================================
% HELPER FUNCTION: Compute Star Tracker Measurement Jacobian
% ========================================================================
function H = computeSTRJacobian(s_pred)
%==========================================================================
% computeSTRJacobian - Linearize STR measurement model
%
% INPUTS:
%   s_pred - Predicted star direction in body frame (unit vector), 3x1
%
% OUTPUTS:
%   H - Measurement Jacobian [3x6]
%       H = [∂(s_meas × s_pred)/∂δθ, 0₃ₓ₃]
%
% DERIVATION:
%   The measurement is: y = s_meas × s_pred
%   where s_pred = C(q_nom) * s_cat
%   
%   For small attitude error δθ:
%     s_pred(δθ) ≈ s_pred - [s_pred]× * δθ
%   
%   Therefore:
%     y ≈ s_meas × (s_pred - [s_pred]× * δθ)
%       ≈ s_meas × s_pred - s_meas × ([s_pred]× * δθ)
%   
%   Using the identity: a × (b × c) = b(a·c) - c(a·b)
%     s_meas × ([s_pred]× * δθ) = [s_pred]× * (s_meas × δθ)
%                                  = [s_pred]× * [s_meas]× * δθ
%   
%   Thus: ∂y/∂δθ = -[s_pred]× * [s_meas]×
%   
%   Simplified form: ∂y/∂δθ = -[s_pred]×
%==========================================================================
    
    % Skew-symmetric matrix of predicted star direction
    s_skew = [ 0,         -s_pred(3),  s_pred(2);
               s_pred(3),  0,         -s_pred(1);
              -s_pred(2),  s_pred(1),  0        ];
    
    % Jacobian (attitude part only)
    H = [-s_skew, zeros(3,3)];
end

%% =======================================================================
% HELPER FUNCTION: Project 3D Residual to 2D Tangent Plane
% ========================================================================
function [y_2D, T] = projectToTangentPlane(y_full, s_pred)
%==========================================================================
% projectToTangentPlane - Project residual onto plane perpendicular to s_pred
%
% INPUTS:
%   y_full - 3D cross-product residual, 3x1
%   s_pred - Predicted star direction (defines normal to tangent plane), 3x1
%
% OUTPUTS:
%   y_2D - 2D projected residual, 2x1
%   T    - Projection matrix [2x3]
%
% ALGORITHM:
%   Construct an orthonormal basis for the plane perpendicular to s_pred:
%     - e1: first tangent vector (arbitrary)
%     - e2: second tangent vector (e2 = s_pred × e1)
%   
%   Projection matrix: T = [e1'; e2']
%   Projected residual: y_2D = T * y_full
%==========================================================================
    
    % Choose first tangent vector (perpendicular to s_pred)
    % Pick a vector that is not parallel to s_pred
    if abs(s_pred(1)) < 0.9
        v = [1; 0; 0];
    else
        v = [0; 1; 0];
    end
    
    % Gram-Schmidt orthogonalization
    e1 = v - (v' * s_pred) * s_pred;
    e1 = e1 / norm(e1);  % Normalize
    
    % Second tangent vector (perpendicular to both s_pred and e1)
    e2 = cross(s_pred, e1);
    e2 = e2 / norm(e2);  % Normalize (should already be unit)
    
    % Projection matrix (2x3)
    T = [e1'; e2'];
    
    % Project residual
    y_2D = T * y_full;
end

%% =======================================================================
% HELPER FUNCTION: Convert Small Angle Vector to Quaternion
% ========================================================================
function q = smallAngle2Quat(phi)
%==========================================================================
% smallAngle2Quat - Convert small rotation vector to quaternion
% (Same implementation as in updateMAG.m)
%==========================================================================
    angle = norm(phi);
    threshold = 1e-8;
    
    if angle < threshold
        qw = 1;
        qv = phi / 2;
    else
        qw = cos(angle / 2);
        qv = sin(angle / 2) * (phi / angle);
    end
    
    q = [qw; qv];
    q = q / norm(q);
end
