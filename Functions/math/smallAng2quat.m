function q = smallAng2quat(phi)
%==========================================================================
% smallAng2quat: Converts small attitude error vector to quaternion using
%                exact Rodrigues formula with first-order approximation for
%                near-zero rotations. Used in MEKF error state reset.
%
% Inputs:
%    phi - Small rotation vector (Gibbs vector) [rad] (3x1)
%          Represents rotation axis scaled by rotation angle
%
% Outputs:
%    q   - Unit quaternion (4x1), scalar-first [qw; qx; qy; qz]
%
% Method:
%    For ||φ|| < threshold: q ≈ [1; φ/2] (first-order)
%    Otherwise: q = [cos(θ/2); sin(θ/2)·(φ/||φ||)] (exact Rodrigues)
%    where θ = ||φ||
%==========================================================================

    % --- 1. Compute Rotation Angle ---
    angle     = norm(phi);
    threshold = 1e-8;
    
    % --- 2. Convert to Quaternion ---
    if angle < threshold
        % First-order approximation to avoid division by zero
        qw = 1;
        qv = phi / 2;
    else
        % Exact Rodrigues formula
        qw = cos(angle / 2);
        qv = sin(angle / 2) * (phi / angle);
    end
    
    % --- 3. Assemble and Normalize ---
    q = [qw; qv];
    q = q / norm(q);

end
