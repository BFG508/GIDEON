function q_inv = quatinv(q)
%==========================================================================
% quatinv: Compute the multiplicative inverse of a quaternion.
%
% Inputs:
%   q - Quaternion in scalar-first convention (4x1):
%       q = [qw; qx; qy; qz].
%
% Outputs:
%   q_inv - Inverse quaternion (4x1) such that:
%           q * q_inv = q_inv * q = [1; 0; 0; 0] (identity quaternion)
%==========================================================================

    % Compute quaternion conjugate: flip sign of vector part
    % q* = [qw; -qx; -qy; -qz]
    q_conj = [q(1); -q(2:4)];
    
    % Compute squared norm: ||q||² = qw² + qx² + qy² + qz²
    % Using norm()^2 is equivalent but more readable than dot(q, q)
    q_norm_squared = norm(q)^2;
    
    % Compute inverse: q_inv = q* / ||q||²
    % For unit quaternions (||q|| = 1), this reduces to q* (conjugate)
    q_inv = q_conj / q_norm_squared;
end
