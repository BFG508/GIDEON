function q = dcm2quat(DCM)
%==========================================================================
% dcm2quat_custom: Convert Direction Cosine Matrix to unit quaternion.
%
% Inputs:
%   DCM - Direction Cosine Matrix (3x3), orthonormal
%
% Outputs:
%   q - Unit quaternion [qw; qx; qy; qz] in scalar-first convention
%
% Method: Shepperd's algorithm (numerically stable)
%==========================================================================

    % Trace of DCM
    trace_DCM = trace(DCM);
    
    % Shepperd's method: select largest component to avoid division by small number
    if trace_DCM > 0
        % qw is largest
        s = sqrt(1 + trace_DCM) * 2;
        qw = 0.25 * s;
        qx = (DCM(3,2) - DCM(2,3)) / s;
        qy = (DCM(1,3) - DCM(3,1)) / s;
        qz = (DCM(2,1) - DCM(1,2)) / s;
    elseif (DCM(1,1) > DCM(2,2)) && (DCM(1,1) > DCM(3,3))
        % qx is largest
        s = sqrt(1 + DCM(1,1) - DCM(2,2) - DCM(3,3)) * 2;
        qw = (DCM(3,2) - DCM(2,3)) / s;
        qx = 0.25 * s;
        qy = (DCM(1,2) + DCM(2,1)) / s;
        qz = (DCM(1,3) + DCM(3,1)) / s;
    elseif (DCM(2,2) > DCM(3,3))
        % qy is largest
        s = sqrt(1 + DCM(2,2) - DCM(1,1) - DCM(3,3)) * 2;
        qw = (DCM(1,3) - DCM(3,1)) / s;
        qx = (DCM(1,2) + DCM(2,1)) / s;
        qy = 0.25 * s;
        qz = (DCM(2,3) + DCM(3,2)) / s;
    else
        % qz is largest
        s = sqrt(1 + DCM(3,3) - DCM(1,1) - DCM(2,2)) * 2;
        qw = (DCM(2,1) - DCM(1,2)) / s;
        qx = (DCM(1,3) + DCM(3,1)) / s;
        qy = (DCM(2,3) + DCM(3,2)) / s;
        qz = 0.25 * s;
    end
    
    % Assemble quaternion
    q = [qw; qx; qy; qz];
    
    % Normalize (should already be unit, but for numerical safety)
    q = q / norm(q);
    
end
