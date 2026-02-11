function Om = skewMatrix(omegaBody)
%==========================================================================
% skewMatrix - Construct Omega matrix for quaternion kinematics
%
% DESCRIPTION:
%   Constructs the skew-symmetric Omega matrix used in the quaternion 
%   differential equation for attitude propagation. This version assumes
%   a scalar-first quaternion convention: q = [qw; qx; qy; qz].
%
% INPUTS:
%   omegaBody - Angular velocity vector in BODY frame [rad/s], 3x1
%
% OUTPUTS:
%   Om        - Omega matrix [4x4] for quaternion kinematics
%
% KINEMATICS EQUATION:
%   The quaternion time derivative is given by:
%     qDot = 1/2 * Omega(omega_b) * q
%
%   where q is the scalar-first quaternion [qw; qx; qy; qz] representing
%   the orientation from inertial to body frame.
%==========================================================================

    wx = omegaBody(1); 
    wy = omegaBody(2); 
    wz = omegaBody(3);
    
    % Construct skew-symmetric Omega matrix (scalar-first convention)
    Om = [ 0, -wx, -wy, -wz;
          wx,   0,  wz, -wy;
          wy, -wz,   0,  wx;
          wz,  wy, -wx,   0];
end
