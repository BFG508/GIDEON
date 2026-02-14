function dy = propragatePerturbedOrbit(t, y, muEarth, REarth, J2Earth, wEarth, inv_BC, cSRP, epoch)
%==========================================================================
% propragatePerturbedOrbit: computes the state derivative for a spacecraft
%                           subject to the following perturbations:
%                             1. J2 Zonal Harmonic (Earth oblateness)
%                             2. Atmospheric Drag (NRLMSISE-00 model)
%                             3. Solar Radiation Pressure (SRP) with cylindrical shadow model
%
% INPUTS:
%   t       - Simulation time [s] from epoch0
%   y       - State vector [6x1]: [rx; ry; rz; vx; vy; vz] (ECI frame)
%   muEarth - Earth's gravitational parameter [m^3/s^2]
%   REarth  - Earth's equatorial radius [m]
%   J2Earth - Second zonal harmonic coefficient [-]
%   wEarth  - Earth's rotation rate [rad/s]
%   inv_BC  - Inverse Ballistic Coefficient for Drag [m^2/kg] * Cd
%   cSRP    - SRP force coefficient [N*m^2/kg] (includes Cr, Area, Mass)
%   epoch   - Simulation start epoch (datetime UTC)
%
% OUTPUTS:
%   dy      - State derivative [6x1]: [vx; vy; vz; ax; ay; az]
%
% ALGORITHM:
%   1. Two-body acceleration (Keplerian)
%   2. J2 Perturbation: Closed-form expression in Cartesian coordinates
%   3. Atmospheric Drag:
%      - Transforms position to LLA
%      - Queries NRLMSISE-00 for density based on solar flux (F10.7)
%      - Computes relative velocity
%      - aDrag = -1/2 * rho * (Cd*A/m) * V_rel^2 * unit(V_rel)
%   4. SRP:
%      - Computes Sun vector
%      - Checks eclipse condition (cylindrical shadow model)
%      - aSRP = -nu * PSun * (Cr*A/m) * unit(Sun-Sat)
%==========================================================================

    rVec = y(1:3);
    vVec = y(4:6);
    rNorm = norm(rVec);
    
    xPos = rVec(1); 
    yPos = rVec(2); 
    zPos = rVec(3);
    
    %% --- 1. J2 Gravity Perturbation ---
    % Acceleration due to central body
    aG = -muEarth * rVec / rNorm^3;
    
    % Acceleration due to J2 (Oblateness)
    % Factor common to all components
    J2Const  = (1.5 * J2Earth * muEarth * REarth^2) / rNorm^5;
    zSquared = (zPos / rNorm)^2;
    
    aJ2   = zeros(3,1);
    aJ2(1) = J2Const * xPos * (5 * zSquared - 1); % x-component
    aJ2(2) = J2Const * yPos * (5 * zSquared - 1); % y-component
    aJ2(3) = J2Const * zPos * (5 * zSquared - 3); % z-component
    
    %% --- 2. Atmospheric Drag (NRLMSISE-00) ---
    % Geodetic altitude check
    % NRLMSISE-00 is valid up to ~ 1000 km, negligible above
    
    % Current simulation time
    currentTime = epoch + seconds(t);
    
    % Transform ECI position to ECEF
    GMST = greenwichSiderealTime(currentTime);
    cosGMST = cos(GMST);
    sinGMST = sin(GMST);
    
    rECEF = [ cosGMST * xPos + sinGMST * yPos;
             -sinGMST * xPos + cosGMST * yPos;
                                         zPos ];
               
    % Convert ECEF to LLA
    LLA = ecef2lla(rECEF', 'WGS84');
    lat = LLA(1);
    lon = LLA(2);
    alt = LLA(3);
    
    if alt > 1000e3
        % Drag is negligible above 1000 km
        aDrag = [0; 0; 0];
    else
        % --- Atmosphere Model Inputs ---
            doy = day(currentTime, 'dayofyear');
        doy_sec = hour(currentTime)*3600 + minute(currentTime)*60 + second(currentTime);
        
        % Solar Flux / Geomagnetic Indices (Static moderate activity)
        % F10.7: Solar Radio Flux
        %    Ap: Geomagnetic planetary index
        F107A = 150; 
        F107  = F107A;
        AP    = 4 * ones(1, 7); % Array of 7 index values
        
        % Obtain density from NRLMSISE-00
        % Output: [He, O, N2, O2, Ar, TotalMassDensity, H, N, AnomalousO]
        [~, rhoOut] = atmosnrlmsise00(alt, lat, lon, year(currentTime), ...
                                      doy, doy_sec, F107A, F107, AP, 'Oxygen');
        rho = rhoOut(6); % Total mass density [kg/m^3]
        
        % Compute Relative Velocity
        vAtm = [-wEarth * yPos; 
                  wEarth * xPos; 
                             0]; 
        vRel     = vVec - vAtm;
        
        % Drag Acceleration
        aDrag = -1/2 * rho * inv_BC * norm(vRel) * vRel;
    end
    
    %% --- 3. Solar Radiation Pressure (SRP) ---
    % Sun Position Vector
    rSun = sunPosition(currentTime);
    
    % Spacecraft-to-Sun Vector
    sVec  = rSun - rVec;
    sDist = norm(sVec);
    sHat  = sVec / sDist; % Unit vector pointing to Sun
    
    % 3. Eclipse Condition (Cylindrical Shadow Model)
    % Project spacecraft position onto Sun direction
    proj = dot(rVec, sHat); 
    
    nu = 1; % Sunlight factor (1 = Sun, 0 = Shadow)
    
    if proj < 0 
        % Spacecraft is on the 'night side' hemisphere
        % Check perpendicular distance to Sun-Earth line
        distPerp = norm(rVec - proj * sHat);
        
        if distPerp < REarth
            % Within Earth's cylindrical shadow
            nu = 0; 
        end
    end
    
    % SRP Acceleration
    if nu > 0
        aSRP = -nu * (cSRP / sDist^2) * sHat;
    else
        aSRP = [0;0;0];
    end
    
    %% --- 4. Total State Derivative ---
    dy = [                   vVec; 
          aG + aJ2 + aDrag + aSRP];
end
