function qNew = quatint(q, omega, dt)
%==========================================================================
% quatint: Integrates quaternion kinematics for attitude propagation using
%          exact exponential map method. Handles small-angle cases to avoid
%          numerical singularities when angular velocity approaches zero.
%
% Inputs:
%    q     - Current quaternion (4x1), scalar-first [qw; qx; qy; qz]
%    omega - Angular velocity vector (body frame) [rad/s] (3x1)
%    dt    - Time step [s]
%
% Outputs:
%    qNew  - Updated quaternion (normalized) (4x1)
%
% Method:
%    Exact exponential map: q(t+dt) = exp(1/2*dt*Ω(ω)) ⊗ q(t)
%    For constant angular velocity: exp(Ω) = [cos(θ/2); sin(θ/2)*ω/||ω||]
%    where θ = ||ω||*dt
%
% REFERENCE:
%    Markley & Crassidis, "Fundamentals of Spacecraft Attitude 
%    Determination and Control", Chapter 6
%==========================================================================

    % --- 1. Normalize input quaternion ---
    q = q / norm(q);
    
    % --- 2. Compute rotation angle ---
    theta     = norm(omega) * dt;
    threshold = 1e-8;
    
    % --- 3. Compute Delta Quaternion ---
    if theta < threshold
        % Small-angle approximation to avoid division by zero
        deltaQ = [1; 1/2 * dt * omega];
    else
        % Exact exponential map
        halfTheta = theta / 2;
        omegaNorm = omega / norm(omega);
        
        deltaQ = [cos(halfTheta);
                  sin(halfTheta) * omegaNorm];
    end
    
    % --- 4. Apply rtation and normalize ---
    qNew = quatmultiply(q, deltaQ);
    qNew = qNew / norm(qNew);

end
