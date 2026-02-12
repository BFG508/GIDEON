function Phi = computeSTM(omega, dt, method)
%==========================================================================
% computeSTM - Compute State Transition Matrix (STM) for MEKF error propagation
%
% INPUTS:
%   omega  - Corrected angular velocity (body frame) [rad/s], 3x1
%   dt     - Time step [s], scalar
%   method - (OPTIONAL) Computation method:
%            'discrete' - First-order discrete approximation (default)
%            'exact'    - Exact closed-form solution (van Loan method)
%
% OUTPUTS:
%   Phi - State Transition Matrix [6x6]
%         Propagates error state: x(k+1) = Phi * x(k) + w(k)
%         where x = [δθ (3x1); b_gyro (3x1)]
%
% MEKF ERROR STATE DYNAMICS (Continuous-Time):
%   dx/dt = F * x + G * w
%   
%   F = [−[ω]×  −I₃]    G = [−I₃   0₃]
%       [ 0₃     0₃ ]        [ 0₃   I₃]
%
%   where [ω]× is the skew-symmetric matrix of angular velocity
%
% DISCRETIZATION:
%   Phi = exp(F * dt) ≈ I + F*dt (first-order)
%   or exact matrix exponential for 'exact' method
%
% REFERENCE:
%   Markley & Crassidis, "Fundamentals of Spacecraft Attitude 
%   Determination and Control", Section 7.5
%==========================================================================

    %% Handle optional method argument
    if nargin < 3
        method = 'discrete';  % Default: first-order approximation
    end
    
    %% Select computation method
    switch lower(method)
        case 'discrete'
            Phi = computeDiscreteSTM(omega, dt);
            
        case 'exact'
            Phi = computeExactSTM(omega, dt);
            
        otherwise
            error('Unknown STM method: %s. Use ''discrete'' or ''exact''.', method);
    end
end

%% =======================================================================
% METHOD 1: First-Order Discrete Approximation
% ========================================================================
function Phi = computeDiscreteSTM(omega, dt)
%==========================================================================
% First-order discrete approximation: Phi ≈ I + F*dt
% Valid for small dt (typically dt < 0.1 s with moderate angular rates)
% Fast and sufficient for most MEKF applications at high sampling rates.
%==========================================================================
    
    % Skew-symmetric matrix of angular velocity
    omega_skew = skewMatrix(omega);
    
    % Continuous-time system matrix F
    % F = [−[ω]×  −I₃]
    %     [ 0₃     0₃ ]
    F = [-omega_skew, -eye(3);
         zeros(3,3),   zeros(3,3)];
    
    % First-order approximation: Phi = I + F*dt
    Phi = eye(6) + F * dt;
end

%% =======================================================================
% METHOD 2: Exact Matrix Exponential (Van Loan Method)
% ========================================================================
function Phi = computeExactSTM(omega, dt)
%==========================================================================
% Exact computation of Phi = exp(F*dt) using closed-form solution.
% More accurate for larger time steps or high angular rates.
%
% For the MEKF, we can derive a closed-form solution:
%   Phi_θθ = I - [ω*dt]× * sin(θ)/θ + [ω*dt]×² * (1-cos(θ))/θ²
%   Phi_θb = -I * dt * sinc(θ/2)
%   Phi_bθ = 0
%   Phi_bb = I
% where θ = ||ω|| * dt
%==========================================================================
    
    % Rotation angle over time step
    theta = norm(omega) * dt;  % [rad]
    
    % Small-angle threshold
    threshold = 1e-8;
    
    if theta < threshold
        %% Small-angle approximation (same as first-order)
        Phi = computeDiscreteSTM(omega, dt);
    else
        %% Exact closed-form solution
        
        % Skew-symmetric matrix of (omega * dt)
        omega_dt_skew = skewMatrix(omega * dt);
        
        % Rodrigues rotation formula for Phi_θθ
        % Phi_θθ = I - [ω*dt]× * sin(θ)/θ + [ω*dt]×² * (1-cos(θ))/θ²
        sinc_theta = sin(theta) / theta;
        one_minus_cos_theta = (1 - cos(theta)) / (theta^2);
        
        Phi_theta_theta = eye(3) - omega_dt_skew * sinc_theta + ...
                          (omega_dt_skew * omega_dt_skew) * one_minus_cos_theta;
        
        % Coupling term Phi_θb (integral of rotation matrix)
        % Phi_θb = -I * dt * sinc(θ/2)
        % For more accuracy, use: -I * dt * (1 - θ²/24 + θ⁴/1920 - ...)
        if theta < 0.1
            % Taylor series for better numerical stability
            sinc_half = 1 - (theta^2)/24 + (theta^4)/1920;
        else
            sinc_half = sin(theta/2) / (theta/2);
        end
        
        Phi_theta_bias = -eye(3) * dt * sinc_half;
        
        % Bias dynamics (random walk, no coupling)
        Phi_bias_theta = zeros(3,3);
        Phi_bias_bias  = eye(3);
        
        % Assemble full STM
        Phi = [Phi_theta_theta, Phi_theta_bias;
               Phi_bias_theta,  Phi_bias_bias];
    end
end

%% =======================================================================
% HELPER FUNCTION: Skew-Symmetric Matrix
% ========================================================================
function S = skewMatrix(v)
%==========================================================================
% skewMatrix - Construct 3x3 skew-symmetric matrix from 3D vector
%
% INPUT:
%   v - 3D vector [v1; v2; v3]
%
% OUTPUT:
%   S - 3x3 skew-symmetric matrix such that S*w = v × w
%
% MATRIX FORM:
%   [v]× = [  0   -v3   v2 ]
%          [ v3    0   -v1 ]
%          [-v2   v1    0  ]
%==========================================================================
    
    S = [ 0,    -v(3),  v(2);
          v(3),  0,    -v(1);
         -v(2),  v(1),  0   ];
end
