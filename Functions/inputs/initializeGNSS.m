function GNSS = initializeGNSS()
%==========================================================================
% initializeGNSS - Defines the hardware and noise parameters for the GNSS
%                  receiver (Space-grade PVT sensor).
%
% OUTPUTS:
%   GNSS - Structure containing sample rate, lever arm, white noise, and 
%          Gauss-Markov parameters for position and velocity correlated errors.
%==========================================================================

    % Update rate [s] (Typically 1 Hz for space receivers)
    GNSS.dt = 1.0; 

    % Lever Arm: Distance from S/C Center of Mass to GNSS Antenna in Body frame
    % Example: Antenna mounted on +Z panel, offset in X
    GNSS.leverArm = [0.5; 0.0; 1.2]; % [m]

    % --- Position Error Model ---
    GNSS.sigmaPosWhite = 1.5;   % High-frequency thermal noise (1-sigma) [m]
    GNSS.tauPos        = 600;   % Correlation time for Ephemeris/Iono errors [s]
    GNSS.sigmaPosGM    = 2.5;   % Correlated noise steady-state (1-sigma) [m]

    % --- Velocity Error Model ---
    GNSS.sigmaVelWhite = 0.05;  % High-frequency noise (1-sigma) [m/s]
    GNSS.tauVel        = 600;   % Correlation time for Clock Wander [s]
    GNSS.sigmaVelGM    = 0.02;  % Correlated noise steady-state (1-sigma) [m/s]

end