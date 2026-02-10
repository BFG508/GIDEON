function euler = dcm2euler_ZYX(DCM)
%==========================================================================
% dcm2euler_ZYX: Extract Euler angles from Direction Cosine Matrix
%                using ZYX (Yaw-Pitch-Roll) convention.
%
% Inputs:
%   DCM   - Direction Cosine Matrix (3x3), representing rotation from 
%           reference frame to body frame. Must be orthonormal (det=1).
%
% Outputs:
%   euler - Euler angles [yaw; pitch; roll] in radians (3x1).
%           Yaw   (ψ): Rotation about Z-axis [-π,   π]
%           Pitch (θ): Rotation about Y-axis [-π/2, π/2]
%           Roll  (φ): Rotation about X-axis [-π,   π]
%
% Gimbal Lock:
%   When |pitch| ≈ 90°, the solution is non-unique (gimbal lock).
%   In this case, roll is set to 0 and only the sum/difference of
%   yaw and roll can be determined.
%==========================================================================

    % Extract pitch angle from DCM(3,1) element
    % Clamp to [-1, 1] to handle numerical errors near ±90°
    pitch = asin(max(min(-DCM(3,1), 1), -1));
    
    % Check for gimbal lock condition (pitch near ±90°)
    if abs(cos(pitch)) > 1e-6
        % Normal case: cos(pitch) is significantly non-zero
        % Extract roll and yaw using standard formulas
        roll = atan2(DCM(3,2), DCM(3,3));
        yaw  = atan2(DCM(2,1), DCM(1,1));
    else
        % Gimbal lock case: pitch ≈ ±90°
        % Yaw and roll are coupled; set roll = 0 by convention
        roll = 0;
        yaw  = atan2(-DCM(1,2), DCM(2,2));
    end
    
    % Return Euler angles as column vector [yaw; pitch; roll]
    euler = [yaw; pitch; roll];
    
end
