function theta = greenwichSiderealTime(dt)
%==========================================================================
% greenwichSiderealTime: Calculates the Greenwich Mean Sidereal Time angle 
%                        for a given UTC datetime.
%
% INPUTS:
%   dt    - Datetime object (UTC) or vector of datetimes.
%
% OUTPUTS:
%   theta - GMST angle [rad], normalized to [0, 2*pi).
%
% ALGORITHM:
%   Uses the IAU-1982 expression for GMST at 0h UT1, approximated here
%   using UTC time (ignoring dUT1 correction < 0.9s).
%   Formula: theta = GMST0 + omega_Earth * fractional_day
%
%==========================================================================

    % Convert datetime to Julian Date (UTC)
    jd = juliandate(dt);
    
    % Calculate Julian Centuries since J2000.0 epoch (T_UT1)
    % (Using UTC as approximation for UT1)
    T = (jd - 2451545.0) / 36525.0;
    
    % --- GMST at 0h UT1 (Mean Sidereal Time) ---
    % Standard IAU formula coefficients (degrees)
    thetaGMST0 = 280.46061837 + ...
                  360.98564736629 * (jd - 2451545.0) + ...
                  0.000387933 * T.^2 - ...
                  T.^3 / 38710000.0;
    
    % --- Add Earth Rotation for Fractional Day ---
    % 360.985647... deg/day is the Earth's mean rotation rate relative to stars
    fracDay = mod(jd, 1.0);
    
    % Total GMST angle
    % Note: In some formulations, the fractional day is included in the T term above.
    % This explicit addition ensures correct handling of the time-of-day component.
    thetaGMST = thetaGMST0 + 360.98564736629 * fracDay;
    
    % Normalize to [0, 2*pi) radians
    theta = deg2rad(mod(thetaGMST, 360));
end
