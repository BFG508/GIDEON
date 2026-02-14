function [scParams, attParams] = initializeSpacecraft()
%==========================================================================
% initializeSpacecraft - Defines physical and attitude control parameters.
%
% Outputs:
%   scParams  - Spacecraft mass, areas, and aerodynamic/radiation coeffs.
%   attParams - Inertia, magnetic dipole, and PD controller gains.
%==========================================================================

    % --- Spacecraft Physical Parameters ---
    scParams.mass     = 100.0; % Spacecraft mass               [kg]
    scParams.areaDrag = 1.0;   % Cross-sectional area for drag [m²]
    scParams.areaSRP  = 1.2;   % Cross-sectional area for SRP  [m²]
    scParams.Cd       = 2.0;   % Atmospheric drag coeff         [-]
    scParams.Cr       = 1.5;   % Solar radiation pressure coeff [-]

    % --- Attitude Control & Dynamics Parameters ---
    % Inertia tensor [kg·m²]
    attParams.I      = diag(5 + 45 * [rand(), rand(), rand()]); 
    
    % Residual magnetic dipole moment [A·m²]
    attParams.dipole = [-0.1 + 0.2*rand(); 
                        -0.1 + 0.2*rand(); 
                        0.05 + 0.25*rand()];
                    
    % PID Controller gains for Nadir pointing
    attParams.Kp     = 0.1   + 1.9*rand();   % Proportional gain [N·m]
    attParams.Ki     = 0.001 + 0.099*rand(); % Integral gain     [N·m/s]
    attParams.Kd     = 0.5   + 4.5*rand();   % Derivative gain   [N·m·s]
    
    % Pointing mode
    attParams.mode   = 'LVLH'; % Local Vertical Local Horizontal
    
end