function Om = skewOmega_custom(omega_b)
%==========================================================================
% skewOmega_custom - Construct Omega matrix for quaternion kinematics
%
% DESCRIPTION:
%   Constructs the skew-symmetric Omega matrix used in the quaternion 
%   differential equation for attitude propagation. This version assumes
%   a scalar-first quaternion convention: q = [qw; qx; qy; qz].
%
% INPUTS:
%   omega_b - Angular velocity vector in BODY frame [rad/s], 3x1
%             [wx; wy; wz]
%
% OUTPUTS:
%   Om      - Omega matrix [4x4] for quaternion kinematics
%
% KINEMATICS EQUATION:
%   The quaternion time derivative is given by:
%     q_dot = 0.5 * Omega(omega_b) * q
%
%   where q is the scalar-first quaternion [qw; qx; qy; qz] representing
%   the orientation from inertial to body frame.
%==========================================================================

    wx = omega_b(1); 
    wy = omega_b(2); 
    wz = omega_b(3);
    
    % Construct skew-symmetric Omega matrix (scalar-first convention)
    Om = [ 0,  -wx, -wy, -wz;
           wx,  0,  wz, -wy;
           wy, -wz,  0,  wx;
           wz,  wy, -wx,  0 ];
end
