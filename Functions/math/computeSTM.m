function Phi = computeSTM(omega, dt)
%==========================================================================
% computeSTM: Computes State Transition Matrix (STM) for MEKF error state
%             propagation using first-order discrete approximation.
%
% Inputs:
%    omega - Corrected angular velocity (body frame) [rad/s] (3x1)
%    dt    - Time step [s]
%
% Outputs:
%    Phi   - State Transition Matrix (6x6)
%            Propagates error state: x(k+1) = Phi*x(k) + w(k)
%            where x = [deltaTheta; biasGyro]
%
% Method:
%    First-order discretization: Phi = I + F*dt
%    
%    F = [−[ω]×  −I₃]
%        [ 0₃     0₃ ]
%    
%    where [ω]× is the skew-symmetric matrix of angular velocity
%
% Reference:
%    Markley & Crassidis, "Fundamentals of Spacecraft Attitude 
%    Determination and Control", Chapter 7
%==========================================================================

    % --- 1. Skew-Symmetric Matrix of Angular Velocity ---
    omegaSkew = [        0, -omega(3),  omega(2);
                  omega(3),         0, -omega(1);
                 -omega(2),  omega(1),        0];
    
    % --- 2. Continuous-Time System Matrix F ---
    F = [-omegaSkew,    -eye(3);
         zeros(3,3), zeros(3,3)];
    
    % --- 3. First-Order Discrete Approximation ---
    Phi = eye(6) + F * dt;

end
