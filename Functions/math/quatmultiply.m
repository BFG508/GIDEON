function qProd = quatmultiply(q1, q2)
%==========================================================================
% quatmultiply: Compute the Hamilton product of two quaternions.
%
% Inputs:
%   q1 - First quaternion in scalar-first convention (4x1):
%        q1 = [qw1; qx1; qy1; qz1] where qw1 is the scalar part and
%        [qx1; qy1; qz1] is the vector part.
%   q2 - Second quaternion in scalar-first convention (4x1):
%        q2 = [qw2; qx2; qy2; qz2] where qw2 is the scalar part and
%        [qx2; qy2; qz2] is the vector part.
%
% Outputs:
%   q_prod - Product quaternion q1 ⊗ q2 in scalar-first convention (4x1).
%            Represents the composition of two rotations: first q2, then q1.
%            WARNING: Quaternion multiplication is non-commutative!
%            q1 ⊗ q2 ≠ q2 ⊗ q1 in general.
%
% Mathematical Definition:
%   Given q1 = (qw1, qv1) and q2 = (qw2, qv2), the Hamilton product is:
%       q_prod = (qw1*qw2 - qv1·qv2,  qw1*qv2 + qw2*qv1 + qv1 × qv2)
%==========================================================================

    % Extract scalar and vector parts from first quaternion
    % CRITICAL: This ordering must match all other quaternion functions
    qw1 = q1(1);    % Scalar part (real component)
    qv1 = q1(2:4);  % Vector part [qx1; qy1; qz1]
    
    % Extract scalar and vector parts from second quaternion
    qw2 = q2(1);    % Scalar part (real component)
    qv2 = q2(2:4);  % Vector part [qx2; qy2; qz2]
    
    % Compute scalar part of product using Hamilton product formula
    % Scalar component: qw1*qw2 - (qv1 · qv2)
    qw = qw1*qw2 - dot(qv1, qv2);
    
    % Compute vector part of product using Hamilton product formula
    % Vector component: qw1*qv2 + qw2*qv1 + (qv1 × qv2)
    qv = qw1*qv2 + qw2*qv1 + cross(qv1, qv2);
    
    % Assemble output quaternion in scalar-first convention
    qProd = [qw; qv];
end
