function q_new = integrateQuaternion(q, omega, dt, method)
%==========================================================================
% integrateQuaternion - Integrate quaternion kinematics for attitude propagation
%
% INPUTS:
%   q      - Current quaternion [qw; qx; qy; qz] (scalar-first), 4x1
%   omega  - Angular velocity vector (body frame) [rad/s], 3x1
%   dt     - Time step [s], scalar
%   method - (OPTIONAL) Integration method:
%            'euler'   - First-order Euler (fast, less accurate)
%            'rk4'     - Fourth-order Runge-Kutta (default, good balance)
%            'exact'   - Exact exponential map (best for large dt/omega)
%
% OUTPUTS:
%   q_new - Updated quaternion (normalized), 4x1
%
% QUATERNION KINEMATICS:
%   dq/dt = 0.5 * Ω(ω) * q
%   where Ω(ω) is the skew-symmetric matrix:
%       Ω(ω) = [  0   -ωx  -ωy  -ωz ]
%              [ ωx    0    ωz  -ωy ]
%              [ ωy   -ωz   0    ωx ]
%              [ ωz    ωy  -ωx   0  ]
%
% REFERENCE:
%   Markley & Crassidis, "Fundamentals of Spacecraft Attitude 
%   Determination and Control", Chapter 7
%==========================================================================

    %% Handle optional method argument
    if nargin < 4
        method = 'rk4';  % Default: 4th-order Runge-Kutta
    end
    
    %% Normalize input quaternion (for numerical safety)
    q = q / norm(q);
    
    %% Select integration method
    switch lower(method)
        case 'euler'
            q_new = integrateEuler(q, omega, dt);
            
        case 'rk4'
            q_new = integrateRK4(q, omega, dt);
            
        case 'exact'
            q_new = integrateExact(q, omega, dt);
            
        otherwise
            error('Unknown integration method: %s. Use ''euler'', ''rk4'', or ''exact''.', method);
    end
    
    %% Normalize output quaternion
    q_new = q_new / norm(q_new);
end

%% =======================================================================
% METHOD 1: First-Order Euler Integration
% ========================================================================
function q_new = integrateEuler(q, omega, dt)
%==========================================================================
% First-order Euler integration: q(t+dt) ≈ q(t) + dt * dq/dt
% Fast but least accurate. Suitable for small dt and low angular rates.
%==========================================================================
    
    % Quaternion derivative
    qDot = 0.5 * omegaMatrix(omega) * q;
    
    % Euler step
    q_new = q + dt * qDot;
end

%% =======================================================================
% METHOD 2: Fourth-Order Runge-Kutta Integration
% ========================================================================
function q_new = integrateRK4(q, omega, dt)
%==========================================================================
% Fourth-order Runge-Kutta (RK4) integration.
% Good balance between accuracy and computational cost.
% Assumes constant angular velocity over [t, t+dt].
%==========================================================================
    
    % k1 = f(t, q)
    k1 = 0.5 * omegaMatrix(omega) * q;
    
    % k2 = f(t + dt/2, q + dt/2 * k1)
    q2 = q + 0.5 * dt * k1;
    q2 = q2 / norm(q2);  % Renormalize intermediate step
    k2 = 0.5 * omegaMatrix(omega) * q2;
    
    % k3 = f(t + dt/2, q + dt/2 * k2)
    q3 = q + 0.5 * dt * k2;
    q3 = q3 / norm(q3);
    k3 = 0.5 * omegaMatrix(omega) * q3;
    
    % k4 = f(t + dt, q + dt * k3)
    q4 = q + dt * k3;
    q4 = q4 / norm(q4);
    k4 = 0.5 * omegaMatrix(omega) * q4;
    
    % Weighted average
    q_new = q + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
end

%% =======================================================================
% METHOD 3: Exact Exponential Map Integration
% ========================================================================
function q_new = integrateExact(q, omega, dt)
%==========================================================================
% Exact solution using quaternion exponential map.
% Most accurate method, especially for large rotation angles.
% Computes: q(t+dt) = exp(0.5 * dt * Ω(ω)) * q(t)
%
% For constant angular velocity:
%   exp(0.5 * dt * Ω(ω)) = [cos(θ/2); sin(θ/2) * ω/||ω||]
%   where θ = ||ω|| * dt
%==========================================================================
    
    % Rotation angle over time step
    theta = norm(omega) * dt;  % [rad]
    
    % Small-angle threshold (below this, use first-order approximation)
    threshold = 1e-8;
    
    if theta < threshold
        % Small-angle approximation: exp(Ω*dt) ≈ I + Ω*dt
        % This avoids division by zero when omega ≈ 0
        delta_q = [1; 0.5 * dt * omega];
    else
        % Exact exponential map
        half_theta = theta / 2;
        omega_norm = omega / norm(omega);  % Rotation axis (unit vector)
        
        delta_q = [cos(half_theta);
                   sin(half_theta) * omega_norm];
    end
    
    % Apply rotation: q_new = delta_q ⊗ q
    q_new = quatmultiply(delta_q, q);
end

%% =======================================================================
% HELPER FUNCTION: Omega Matrix (Quaternion Kinematics)
% ========================================================================
function Omega = omegaMatrix(omega)
%==========================================================================
% omegaMatrix - Construct skew-symmetric matrix for quaternion kinematics
%
% INPUT:
%   omega - Angular velocity [rad/s], 3x1
%
% OUTPUT:
%   Omega - 4x4 skew-symmetric matrix such that dq/dt = 0.5 * Omega * q
%
% MATRIX FORM:
%   Ω(ω) = [  0   -ωx  -ωy  -ωz ]
%          [ ωx    0    ωz  -ωy ]
%          [ ωy   -ωz   0    ωx ]
%          [ ωz    ωy  -ωx   0  ]
%==========================================================================
    
    wx = omega(1);
    wy = omega(2);
    wz = omega(3);
    
    Omega = [ 0,  -wx, -wy, -wz;
             wx,   0,   wz, -wy;
             wy,  -wz,  0,   wx;
             wz,   wy, -wx,   0 ];
end

%% =======================================================================
% HELPER FUNCTION: Quaternion Multiplication
% ========================================================================
function q = quatmultiply(q1, q2)
%==========================================================================
% quatmultiply - Multiply two quaternions (scalar-first convention)
%
% INPUTS:
%   q1, q2 - Quaternions [qw; qx; qy; qz], 4x1
%
% OUTPUT:
%   q - Product quaternion q1 ⊗ q2, 4x1
%
% FORMULA:
%   q = [q1w*q2w - q1v·q2v;
%        q1w*q2v + q2w*q1v + q1v × q2v]
%==========================================================================
    
    % Extract scalar and vector parts
    q1w = q1(1);
    q1v = q1(2:4);
    
    q2w = q2(1);
    q2v = q2(2:4);
    
    % Quaternion multiplication
    qw = q1w * q2w - dot(q1v, q2v);
    qv = q1w * q2v + q2w * q1v + cross(q1v, q2v);
    
    q = [qw; qv];
end
