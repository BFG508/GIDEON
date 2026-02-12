function [rECI, vECI] = simulateOrbit(t, orbitalElements, scParams)
%==========================================================================
% simulateOrbit: Propagates a spacecraft orbit in the ECI frame using numerical integration
%                (ODE113). The force model includes Earth's J2 oblateness (zonal harmonic),
%                atmospheric drag with exponential density, and solar radiation pressure (SRP).
%                It handles eclipses for SRP and Earth rotation effects for drag.
%
% INPUTS:
%    t               - Time vector [s] (relative to epoch, starting at 0)
%    orbitalElements - Struct with fields:
%        .SMA        - Semi-major axis [m]
%        .ECC        - Eccentricity [0, 1)
%        .INC        - Inclination [deg]
%        .RAAN       - RAAN [deg]
%        .AOP        - Argument of Perigee [deg]
%        .TA         - True Anomaly at t=0 [deg]
%        .epoch      - Datetime object (UTC) for Sun position calculation
%
%    scParams        - Struct with spacecraft physical parameters:
%        .Cd         - Drag coefficient                [-]
%        .Cr         - Radiation pressure coefficient  [-]
%        .areaDrag   - Cross-sectional area for drag  [m²]
%        .areaSRP    - Cross-sectional area for SRP   [m²]
%        .mass       - Spacecraft mass                [kg]
%
% OUTPUTS:
%   rECI           - Position vectors in ECI   [m], 3xN
%   vECI           - Velocity vectors in ECI [m/s], 3xN
%
%==========================================================================

    %% --- 1. Physical Constants (WGS84) ---
    muEarth = 3.986004418e14;  % [m^3/s^2]
    REarth  = 6378137.0;       % [m]
    J2Earth = 1.08262668e-3;   % J2 Zonal Harmonic
    wEarth  = 7.2921150e-5;    % Earth rotation rate [rad/s]
    PSun    = 4.56e-6;         % Solar Radiation Pressure at 1 AU [N/m^2]
    AU      = 149597870700;    % Astronomical Unit [m]
    
    %% --- 2. Initial State Vector (COE -> ECI) ---
    SMA  = orbitalElements.SMA;
    ECC  = orbitalElements.ECC;
    INC  = deg2rad(orbitalElements.INC);
    RAAN = deg2rad(orbitalElements.RAAN);
    AOP  = deg2rad(orbitalElements.AOP);
    TA   = deg2rad(orbitalElements.TA);
    
    % Epoch for sun position and atmospheric model
    epoch0 = orbitalElements.epoch;
    
    p    = SMA * (1 - ECC^2);
    rPQW = (p / (1 + ECC*cos(TA))) * [ cos(TA);         sin(TA); 0];
    vPQW =         sqrt(muEarth/p) * [-sin(TA); (ECC + cos(TA)); 0];
    
    % Rotation Matrix PQW -> ECI
    R_AOP  = [ cos(AOP),  -sin(AOP),         0; 
               sin(AOP),   cos(AOP),         0;
                      0,          0,         1];
    R_INC  = [        1,          0,         0; 
                      0,   cos(INC), -sin(INC); 
                      0,   sin(INC),  cos(INC)];
    R_RAAN = [cos(RAAN), -sin(RAAN),         0; 
              sin(RAAN),  cos(RAAN),         0; 
                      0,          0,         1];
    Q = R_RAAN * R_INC * R_AOP;
    
    r0 = Q * rPQW;
    v0 = Q * vPQW;
    
    %% --- 3. Numerical Integration (ODE113) ---
    y0 = [r0; v0]; 
    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-5);
    
    if numel(t) < 2, rECI = r0'; vECI = v0'; 
        return; 
    end
    
    % Pre-compute ballistic and SRP coefficients
    inv_BC_drag = (scParams.Cd * scParams.areaDrag) / scParams.mass;
    coeff_SRP   = (scParams.Cr *  scParams.areaSRP) / scParams.mass * PSun * AU^2;
    
    % Integrate
    [~, Yout] = ode113(@(tSim, y) propragatePerturbedOrbit(tSim, y, ...
        muEarth, REarth, J2Earth, wEarth, inv_BC_drag, coeff_SRP, epoch0), ...
        t, y0, opts);
    
    rECI = Yout(:, 1:3)';
    vECI = Yout(:, 4:6)';
end