function [r_ECI, v_ECI] = simulateOrbit(t, orbitalElements, scParams)
%==========================================================================
% simulateOrbit - Orbit Propagator: J2 + Drag + Solar Radiation Pressure
%==========================================================================
%
% INPUTS:
%   t               - Time vector [s] (relative to epoch, starting at 0)
%   orbitalElements - Struct with fields:
%       .sma        - Semi-major axis [m]
%       .ecc        - Eccentricity [0, 1)
%       .inc_deg    - Inclination [deg]
%       .raan_deg   - RAAN [deg]
%       .argp_deg   - Argument of Perigee [deg]
%       .nu_deg     - True Anomaly at t=0 [deg]
%       .epoch      - datetime object (UTC) for Sun position calculation
%
%   scParams        - Struct with spacecraft physical parameters:
%       .Cd         - Drag coefficient (typ. 2.2)
%       .Cr         - Radiation pressure coefficient (typ. 1.0-1.8)
%       .Area_drag  - Cross-sectional area for drag [m^2]
%       .Area_srp   - Cross-sectional area for SRP [m^2]
%       .Mass       - Spacecraft mass [kg]
%
% OUTPUTS:
%   r_ECI           - Position vectors in ECI [m], 3xN
%   v_ECI           - Velocity vectors in ECI [m/s], 3xN
%
%==========================================================================

    % --- 1. Physical Constants (WGS84) ---
    mu_Earth = 3.986004418e14;  % [m^3/s^2]
    R_Earth  = 6378137.0;       % [m]
    J2_Earth = 1.08262668e-3;   % J2 Zonal Harmonic
    w_Earth  = 7.2921150e-5;    % Earth rotation rate [rad/s]
    P_sun    = 4.56e-6;         % Solar Radiation Pressure at 1 AU [N/m^2]
    AU       = 149597870700;    % Astronomical Unit [m]
    
    % --- 2. Initial State Vector (COE -> ECI) ---
    SMA  = orbitalElements.sma;
    ECC  = orbitalElements.ecc;
    INC  = deg2rad(orbitalElements.inc_deg);
    RAAN = deg2rad(orbitalElements.raan_deg);
    AOP  = deg2rad(orbitalElements.argp_deg);
    TA   = deg2rad(orbitalElements.nu_deg);
    
    % Epoch for sun position and atmospheric model
    epoch0 = orbitalElements.epoch;
    
    p = SMA * (1 - ECC^2);
    r_pqw = (p / (1 + ECC*cos(TA))) * [cos(TA); sin(TA); 0];
    v_pqw = sqrt(mu_Earth/p) * [-sin(TA); (ECC + cos(TA)); 0];
    
    % Rotation Matrix PQW -> ECI (3-1-3 transposed)
    R_w  = [cos(AOP) -sin(AOP) 0; sin(AOP) cos(AOP) 0; 0 0 1];
    R_i  = [1 0 0; 0 cos(INC) -sin(INC); 0 sin(INC) cos(INC)];
    R_Om = [cos(RAAN) -sin(RAAN) 0; sin(RAAN) cos(RAAN) 0; 0 0 1];
    Q = R_Om * R_i * R_w;
    
    r0 = Q * r_pqw;
    v0 = Q * v_pqw;
    
    % --- 3. Numerical Integration (ODE45) ---
    y0 = [r0; v0]; 
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
    
    if numel(t) < 2, r_ECI = r0'; v_ECI = v0'; return; end
    
    % Pre-compute ballistic/SRP coefficients
    inv_BC_drag = (scParams.Cd * scParams.Area_drag) / scParams.Mass;
    coeff_SRP   = (scParams.Cr * scParams.Area_srp) / scParams.Mass * P_sun * AU^2;
    
    % Integrate
    [~, Y_out] = ode45(@(t_sim, y) dynamicsFull(t_sim, y, ...
        mu_Earth, R_Earth, J2_Earth, w_Earth, inv_BC_drag, coeff_SRP, epoch0), ...
        t, y0, options);
    
    r_ECI = Y_out(:, 1:3)';
    v_ECI = Y_out(:, 4:6)';
end

% --- Internal Dynamics: J2 + Drag + SRP ---
function dy = dynamicsFull(t_sim, y, mu, Re, J2, w_E, inv_BC_drag, coeff_SRP, epoch0)
    r_vec = y(1:3);
    v_vec = y(4:6);
    r_norm = norm(r_vec);
    
    x = r_vec(1); y_pos = r_vec(2); z = r_vec(3);
    
    % --- A. J2 Gravity ---
    a_2body = -mu * r_vec / r_norm^3;
    
    pre = (1.5 * J2 * mu * Re^2) / r_norm^5;
    z_sq = 5 * (z / r_norm)^2;
    
    a_j2 = pre * [ x * (z_sq - 1);
                   y_pos * (z_sq - 1);
                   z * (z_sq - 3) ];
    
    % --- B. Atmospheric Drag ---
    h_m = r_norm - Re;
    
    % Altitude threshold for negligible density (avoid NRLMSISE-00 warnings)
    if h_m > 1000e3
        a_drag = [0;0;0];
    else
        % Convert ECI to geodetic coordinates for NRLMSISE-00
        current_time = epoch0 + seconds(t_sim);
        
        % ECI to ECEF (simple rotation by Greenwich Hour Angle)
        theta_GMST = greenwichSiderealTime(current_time);
        DCM_ECI2ECEF = [cos(theta_GMST)  sin(theta_GMST) 0;
                       -sin(theta_GMST)  cos(theta_GMST) 0;
                        0                0               1];
        r_ECEF = DCM_ECI2ECEF * r_vec;
        
        % ECEF to LLA
        lla = ecef2lla(r_ECEF');
        lat_deg = lla(1);
        lon_deg = lla(2);
        alt_m   = lla(3);
        
        % Only call NRLMSISE-00 for valid altitude range
        if alt_m > 1000e3
            % Above nominal model limit, assume negligible density
            rho = 0;
        else
            % Extract time components for NRLMSISE-00
            yr  = year(current_time);
            doy = day(current_time, 'dayofyear');
            ut_sec = hour(current_time)*3600 + minute(current_time)*60 + second(current_time);
            
            % Solar and geomagnetic activity indices (nominal values for moderate activity)
            F107A = 150;  % 81-day average of F10.7 solar flux [10^-22 W/(m^2·Hz)]
            F107  = 150;  % Daily F10.7 solar flux [10^-22 W/(m^2·Hz)]
            APH   = [4, 4, 4, 4, 4, 4, 4]; % Geomagnetic Ap index (7 values: daily + 3-hourly)
            
            % Call NRLMSISE-00 atmospheric model with activity indices
            [~, rho_out] = atmosnrlmsise00(alt_m, lat_deg, lon_deg, yr, doy, ut_sec, ...
                                            F107A, F107, APH);
            rho = rho_out(6); % Column 6 = total mass density [kg/m³]
        end
        
        % Drag acceleration
        if rho > 0
            % V_rel = V_inertial - (w_earth x r)
            v_atm = [-w_E * y_pos; w_E * x; 0];
            v_rel = v_vec - v_atm;
            v_rel_norm = norm(v_rel);
            
            a_drag = -0.5 * rho * inv_BC_drag * v_rel_norm * v_rel;
        else
            a_drag = [0;0;0];
        end
    end
    
    % --- C. Solar Radiation Pressure (SRP) ---
    % 1. Sun Position at current t_sim
    r_sun = sunPositionApprox(epoch0 + seconds(t_sim));
    
    % 2. Vector from Spacecraft to Sun
    s_vec = r_sun - r_vec;
    s_dist = norm(s_vec);
    s_hat  = s_vec / s_dist;
    
    % 3. Shadow Function (Cylindrical Model)
    proj = dot(r_vec, s_hat); 
    
    nu = 1; % Assume sunlight
    if proj < 0 
        % Spacecraft is 'behind' Earth (relative to Sun)
        dist_perp = norm(r_vec - proj * s_hat);
        if dist_perp < Re
            nu = 0; % In shadow (Eclipse)
        end
    end
    
    % 4. Acceleration
    if nu > 0
        a_srp = -nu * (coeff_SRP / s_dist^2) * s_hat;
    else
        a_srp = [0;0;0];
    end
    
    % --- Total ---
    dy = [v_vec; a_2body + a_j2 + a_drag + a_srp];
end

% --- Helper: Simple Solar Ephemeris ---
function r_sun = sunPositionApprox(dt)
    % Analytical approx for Sun vector in ECI
    % Ref: Montenbruck & Gill, Eq 3.3
    jd = juliandate(dt);
    T = (jd - 2451545.0) / 36525.0;
    
    mean_long = 280.460 + 36000.77*T;
    mean_anom = 357.5277233 + 35999.05*T;
    
    lambda_sun = deg2rad(mean_long + 1.914666471*sin(deg2rad(mean_anom)) + ...
                 0.019994643*sin(deg2rad(2*mean_anom)));
    epsilon    = deg2rad(23.439291 - 0.0130042*T);
    
    R_AU = 1.000140612 - 0.016708617*cos(deg2rad(mean_anom)) - ...
           0.000139589*cos(deg2rad(2*mean_anom));
       
    AU_m = 149597870700;
    r_mag = R_AU * AU_m;
    
    r_sun = r_mag * [cos(lambda_sun); 
                     cos(epsilon)*sin(lambda_sun); 
                     sin(epsilon)*sin(lambda_sun)];
end

% --- Helper: Greenwich Mean Sidereal Time ---
function theta = greenwichSiderealTime(dt)
    % Compute GMST angle from UTC datetime
    % Simplified formula (sufficient for drag calculations)
    jd = juliandate(dt);
    T = (jd - 2451545.0) / 36525.0;
    
    % GMST at 0h UT (degrees)
    theta_GMST0 = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + ...
                  0.000387933 * T^2 - T^3 / 38710000.0;
    
    % Add fraction of day
    frac_day = mod(jd, 1.0);
    theta_GMST = theta_GMST0 + 360.98564736629 * frac_day;
    
    % Normalize to [0, 360) and convert to radians
    theta = deg2rad(mod(theta_GMST, 360));
end
