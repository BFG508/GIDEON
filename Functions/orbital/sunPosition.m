function rSun = sunPosition(dt)
%==========================================================================
% sunPosition: calculates the Sun's position vector (ECI) using a 
%              low-precision analytical series approximation.
%
% INPUTS:
%   dt    - Datetime object (UTC).
%
% OUTPUTS:
%   rSun - Sun position vector in ECI frame [m] (3xN).
%
% ALGORITHM:
%   1. Computes mean longitude and mean anomaly from Julian Date.
%   2. Solves for ecliptic longitude and distance.
%   3. Transforms from Ecliptic to ECI using the obliquity of the ecliptic.
%==========================================================================

    % Convert datetime to Julian Date (UTC)
    jd = juliandate(dt);
    
    % Julian Centuries since J2000.0
    T = (jd - 2451545.0) / 36525.0;
    
    % --- Mean Elements of Earth's Orbit ---
    % Mean Longitude of the Sun [deg]
    sunMeanLon = 280.460 + 36000.77 * T;
    
    % Mean Anomaly of the Sun [deg]
    sunMA     = 357.5277233 + 35999.05 * T;
    
    % --- Ecliptic Coordinates ---
    % Ecliptic Longitude (lambda) via equation of center approximation [rad]
    sunLambda = deg2rad(sunMeanLon + ...
                        1.914666471 * sin(deg2rad(sunMA)) + ...
                        0.019994643 * sin(deg2rad(2 * sunMA)));
    
    % Obliquity of the Ecliptic (epsilon) - Earth's axial tilt [rad]
    epsilon = deg2rad(23.439291 - 0.0130042 * T);
    
    % Distance from Earth to Sun [AU]
    R_AU = 1.000140612 - ...
           0.016708617 * cos(deg2rad(sunMA)) - ...
           0.000139589 * cos(deg2rad(2 * sunMA));
       
    % Astronomical Unit in meters (IAU 2012 constant)
    AU = 149597870700; 
    
    % Magnitude of Sun vector [m]
    rMag = R_AU * AU;
    
    % --- Transform to ECI (Equatorial) Frame ---
    % Rotation from Ecliptic to Equatorial by obliquity epsilon
    % r_x = r * cos(lambda)
    % r_y = r * cos(eps) * sin(lambda)
    % r_z = r * sin(eps) * sin(lambda)
    rSun = rMag .* [cos(sunLambda); 
                      cos(epsilon) .* sin(sunLambda); 
                      sin(epsilon) .* sin(sunLambda)];
end
