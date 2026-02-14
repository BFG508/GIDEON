function B_ECI = simulateGeomagneticField(t, rECI, epoch)
%==========================================================================
% simulateGeomagneticField: Computes the Earth's magnetic field in the ECI 
%                          frame using the IGRF-13 model.
%
% Inputs:
%   t      - Time vector                               [s], 1xN
%   rECI   - Position in ECI frame                     [m], 3xN
%   epoch  - Simulation start UTC time   (datetime object)
%
% Outputs:
%   B_ECI  - Magnetic field in ECI frame              [nT], 3xN
%
% Algorithm:
%   1. Downsamples time to compute ECI-to-ECEF DCM efficiently.
%   2. Rotates ECI positions to ECEF and converts to LLA.
%   3. Queries IGRF-13 for the magnetic field in NED frame using the 
%      exact decimal year derived from the epoch.
%   4. Rotates the magnetic field from NED -> ECEF -> ECI.
%==========================================================================

    N = numel(t);
    
    % 1. Downsample DCM calculation (ECI → ECEF transformation is slow)
    dt_original = mean(diff(t));
    dt_DCM      = 1.0;  % Calculate DCM every 1 second
    downsample_factor = max(1, round(dt_DCM / dt_original));
    
    t_DCM = t(1:downsample_factor:end);
    % Ensure the last point is included to prevent extrapolation issues
    if t_DCM(end) ~= t(end)
        t_DCM = [t_DCM, t(end)]; 
    end
    
    % Compute DCM_ECI2ECEF at downsampled rate
    DCM_ECI2ECEF_ds = dcmeci2ecef('IAU-2000/2006', epoch + seconds(t_DCM'));
    
    % Interpolate DCM to original high-frequency time vector
    DCM_ECI2ECEF = zeros(3, 3, N);
    for i = 1:3
        for j = 1:3
            dcm_ij = squeeze(DCM_ECI2ECEF_ds(i,j,:));
            DCM_ECI2ECEF(i,j,:) = interp1(t_DCM, dcm_ij, t, 'linear', 'extrap');
        end
    end
    
    % 2. Convert ECI to ECEF positions
    rECEF = zeros(3, N);
    for k = 1:N
        rECEF(:,k) = squeeze(DCM_ECI2ECEF(:,:,k)) * rECI(:,k);
    end
    
    % Convert ECEF to Geodetic (LLA) using WGS84
    LLA = ecef2lla(rECEF', 'WGS84');
    lat = LLA(:,1);          % [deg]
    lon = LLA(:,2);          % [deg]
    alt_km = LLA(:,3) / 1e3; %  [km]
    
    % 3. Evaluate IGRF-13 model
    % Compute the exact decimal year from the provided datetime epoch
    yr = year(epoch);
    
    % Number of days in the current year (365 or 366 for leap years)
    daysInYear = day(datetime(yr, 12, 31), 'dayofyear'); 
    
    % Precise decimal year (incorporating days and hours)
    decYear = yr + (day(epoch, 'dayofyear') - 1 + hour(epoch)/24) / daysInYear;
    if decYear > 2025
        decYear = 2025;
    end
    
    XYZ_NED = igrfmagm(alt_km, lat, lon, decYear * ones(N,1), 13);
    
    % 4. Convert NED to ECEF to ECI
    [B_ECEFx, B_ECEFy, B_ECEFz] = ned2ecefv(XYZ_NED(:,1), XYZ_NED(:,2), XYZ_NED(:,3), lat, lon);
    B_ECEF = [B_ECEFx'; B_ECEFy'; B_ECEFz'];
    
    B_ECI = zeros(3, N);
    for k = 1:N
        % Transpose to get ECEF to ECI
        DCM_ECEF2ECI = squeeze(DCM_ECI2ECEF(:,:,k))';
        B_ECI(:,k) = DCM_ECEF2ECI * B_ECEF(:,k);
    end

end