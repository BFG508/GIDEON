function DCM = quat2dcm(q)
%==========================================================================
% quat2dcm: Convert unit quaternion to Direction Cosine Matrix.
%
% Inputs:
%   q - Unit quaternion in scalar-first convention (4x1):
%       q = [qw; qx; qy; qz] = [q0; q1; q2; q3]
%       where qw is the scalar part and [qx; qy; qz] is the vector part.
%       Must satisfy: ||q|| = sqrt(qw^2 + qx^2 + qy^2 + qz^2) = 1
%
% Outputs:
%   DCM - Direction Cosine Matrix (3x3), orthonormal with det(DCM) = 1.
%         Represents passive rotation (coordinate transformation) from
%         reference frame to body frame: vBody = DCM * vReference
%==========================================================================

    % Extract quaternion components using scalar-first convention
    % CRITICAL: This ordering must match all other quaternion functions
    qw = q(1); % Scalar part (real component)
    qx = q(2); % Vector part: i-component
    qy = q(3); % Vector part: j-component
    qz = q(4); % Vector part: k-component
    
    % Construct Direction Cosine Matrix using closed-form quaternion formula
    % Each element is derived from the quaternion product q·v·q*
    % where v is expressed as a pure quaternion [0; v]
    DCM = [1-2*(qy^2  + qz^2),    2*(qx*qy + qz*qw),   2*(qx*qz - qy*qw);
             2*(qx*qy - qz*qw), 1-2*(qx^2  + qz^2),    2*(qy*qz + qx*qw);
             2*(qx*qz + qy*qw),   2*(qy*qz - qx*qw), 1-2*(qx^2  + qy^2)];
end