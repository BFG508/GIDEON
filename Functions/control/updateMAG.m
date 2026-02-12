function ekf = updateMAG(ekf, B_meas, B_ECI)
%==========================================================================
% updateMAG - Magnetometer measurement update for MEKF
%
% INPUTS:
%   ekf    - EKF state structure (from initializeMEKF or predictMEKF)
%   B_meas - Measured magnetic field in body frame [nT], 3x1
%   B_ECI  - True magnetic field in ECI frame [nT], 3x1
%
% OUTPUTS:
%   ekf - Updated EKF state structure with corrected attitude and covariance
%
% ALGORITHM:
%   1. Predict measurement using nominal attitude
%   2. Compute innovation (measurement residual)
%   3. Linearize measurement model (Jacobian H)
%   4. Compute Kalman gain
%   5. Update error state
%   6. Update covariance (Joseph form for numerical stability)
%   7. Apply correction to nominal quaternion and reset error state
%==========================================================================

    %% 1. PREDICT MEASUREMENT
    % Rotate predicted field from ECI to body using nominal attitude
    DCM_nom = quat2dcm(ekf.q_nom);
    B_pred  = DCM_nom * B_ECI;  % [nT]
    
    %% 2. INNOVATION (Measurement Residual)
    y = B_meas - B_pred;  % [nT], 3x1
    
    %% 3. MEASUREMENT JACOBIAN
    % H = [∂h/∂δθ, ∂h/∂b_gyro]
    % For magnetometer: h(x) = C(q_nom) * B_ECI
    % The Jacobian w.r.t. attitude error δθ is: -[B_pred]×
    % The Jacobian w.r.t. gyro bias is zero (bias doesn't affect MAG)
    H = computeMagJacobian(B_pred);  % 3x6
    
    %% 4. KALMAN GAIN
    S = H * ekf.P * H' + ekf.R_mag;  % Innovation covariance
    K = ekf.P * H' / S;              % Kalman gain, 6x3
    
    %% 5. UPDATE ERROR STATE
    ekf.x = ekf.x + K * y;  % a posteriori error state
    
    %% 6. UPDATE COVARIANCE (Joseph Form)
    % Joseph form guarantees positive semi-definiteness even with numerical errors
    I   = eye(6);
    IKH = I - K * H;
    ekf.P = IKH * ekf.P * IKH' + K * ekf.R_mag * K';
    
    %% 7. RESET: Apply Correction to Nominal Quaternion
    % Convert small-angle error to quaternion correction
    delta_q = smallAngle2Quat(ekf.x(1:3));
    
    % Update nominal quaternion: q_nom = δq ⊗ q_nom
    ekf.q_nom = quatmultiply(delta_q, ekf.q_nom);
    ekf.q_nom = ekf.q_nom / norm(ekf.q_nom);  % Renormalize
    
    % Reset attitude error (MEKF standard practice)
    ekf.x(1:3) = zeros(3,1);
    
    % Note: Bias error ekf.x(4:6) is NOT reset
end

%% =======================================================================
% HELPER FUNCTION: Compute Magnetometer Measurement Jacobian
% ========================================================================
function H = computeMagJacobian(B_pred)
%==========================================================================
% computeMagJacobian - Linearize magnetometer measurement model
%
% INPUTS:
%   B_pred - Predicted magnetic field in body frame [nT], 3x1
%
% OUTPUTS:
%   H - Measurement Jacobian [3x6]
%       H = [∂B_body/∂δθ, ∂B_body/∂b_gyro]
%
% DERIVATION:
%   The measurement model is: B_body = C(q) * B_ECI
%   For small attitude error δθ:
%     C(δq ⊗ q_nom) ≈ C(q_nom) * (I - [δθ]×)
%   Therefore:
%     B_body ≈ C(q_nom) * (I - [δθ]×) * C(q_nom)' * B_body_nom
%            ≈ B_pred - [B_pred]× * δθ
%   Thus:
%     ∂B_body/∂δθ = -[B_pred]×
%     ∂B_body/∂b_gyro = 0₃ₓ₃ (bias doesn't affect MAG measurement)
%==========================================================================
    
    % Skew-symmetric matrix of predicted field
    B_skew = [ 0,         -B_pred(3),  B_pred(2);
               B_pred(3),  0,         -B_pred(1);
              -B_pred(2),  B_pred(1),  0        ];
    
    % Jacobian: H = [-[B_pred]×, 0]
    H = [-B_skew, zeros(3,3)];
end

%% =======================================================================
% HELPER FUNCTION: Convert Small Angle Vector to Quaternion
% ========================================================================
function q = smallAngle2Quat(phi)
%==========================================================================
% smallAngle2Quat - Convert small rotation vector to quaternion
%
% INPUTS:
%   phi - Small rotation vector (Gibbs vector or MRP) [rad], 3x1
%         Represents rotation axis scaled by rotation angle
%
% OUTPUTS:
%   q - Unit quaternion [qw; qx; qy; qz] (scalar-first), 4x1
%
% ALGORITHM:
%   For small angles ||phi|| << 1:
%     q ≈ [1; phi/2]  (first-order approximation)
%   
%   For better accuracy (second-order):
%     angle = ||phi||
%     axis  = phi / ||phi||
%     q = [cos(angle/2); sin(angle/2) * axis]
%
% NOTE: We use second-order to avoid singularities when error grows.
%==========================================================================
    
    % Magnitude of rotation
    angle = norm(phi);  % [rad]
    
    % Small-angle threshold (below this, use first-order approx)
    threshold = 1e-8;
    
    if angle < threshold
        % First-order approximation: q ≈ [1; phi/2]
        qw = 1;
        qv = phi / 2;
    else
        % Exact conversion (Rodrigues formula)
        qw = cos(angle / 2);
        qv = sin(angle / 2) * (phi / angle);  % axis * sin(angle/2)
    end
    
    % Assemble quaternion
    q = [qw; qv];
    
    % Normalize (should already be unit, but for numerical safety)
    q = q / norm(q);
end
