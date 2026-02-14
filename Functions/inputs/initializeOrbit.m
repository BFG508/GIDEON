function orbitalElems = initializeOrbit(epoch)
%==========================================================================
% initializeOrbit - Defines the orbital elements for the simulation.
%
% Inputs:
%   epoch           - Datetime object (UTC) representing simulation start
%
% Outputs:
%   orbitalElems    - Struct containing Keplerian elements and epoch
%==========================================================================

    % Default epoch if none is provided
    if nargin < 1
        epoch = datetime(2026, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
    end

    % Orbital elements (randomized Low Earth Orbit - LEO)
    orbitalElems.SMA   = 6378e3 + 200e3 + 400e3 * rand(); % Semi-major axis   [m]
    orbitalElems.ECC   =                  0.01 * rand();  % Eccentricity      [-]
    orbitalElems.INC   =                   180 * rand();  % Inclination     [deg]
    orbitalElems.RAAN  =                   360 * rand();  % Right Ascension [deg]
    orbitalElems.AOP   =                   360 * rand();  % Arg. of Perigee [deg]
    orbitalElems.TA    =                   360 * rand();  % True Anomaly    [deg]
    orbitalElems.epoch = epoch;                           % UTC Epoch
    
end