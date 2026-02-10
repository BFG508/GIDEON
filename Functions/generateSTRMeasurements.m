function meas = generateSTRMeasurements(STR, N_STR, catalog, DCM_true, omega_body)
%==========================================================================
% generateSTRMeasurements: Simulate star tracker observations with
%                          realistic noise models.
%
% Inputs:
%   STR         - Structure array (1xN_STR) with STR specifications:
%                 .DCM_B2S           - Body-to-Sensor DCM (3x3)
%                 .FOV_deg           - Field of view [deg]
%                 .resolution        - Sensor dimensions [pixels]
%                 .pixel_size        - Pixel pitch [m]
%                 .focal_length      - Optical focal length [m]
%                 .centroid_accuracy - Centroiding std dev [pixels] for 
%                                      reference magnitude (magnitude_limit)
%                 .magnitude_limit   - Reference magnitude for centroid_accuracy
%                 .update_rate       - Measurement frequency [Hz]
%   N_STR       - Number of active star trackers (1 or 2).
%   catalog     - Star catalog structure:
%                 .r_ECI             - 3xN unit vectors in ECI frame
%                 .magnitude         - Nx1 absolute magnitudes
%   DCM_true    - True ECI-to-Body DCM (3x3).
%   omega_body  - Satellite angular velocity in body frame [rad/s].
%
% Outputs:
%   meas        - Structure array (1xN_STR) with measurements:
%                 .N_stars           - Number of detected stars
%                 .star_id           - Catalog indices of detected stars
%                 .r_body            - 3xM unit vectors in body frame
%                 .r_ECI_ref         - 3xM reference unit vectors (ECI)
%                 .magnitude         - Mx1 absolute magnitudes
%                 .weights           - Mx1 measurement weights (flux-based)
%                 .pixel_x_noisy     - Mx1 pixel coordinates (for viz)
%                 .pixel_y_noisy     - Mx1 pixel coordinates (for viz)
%
% Method:
%   1. Transform catalog stars to each sensor frame via DCM_B2S * DCM_true.
%   2. Apply gnomonic projection to focal plane (pinhole camera model).
%   3. Filter stars within FOV and physical sensor boundaries.
%   4. Add magnitude-dependent Gaussian centroiding noise and motion blur.
%   5. Back-project noisy pixel coordinates to unit vectors.
%   6. Compute flux-based weights (brighter stars = higher weight).
%
% Notes:
%   - Minimum 3 stars required per STR for attitude determination.
%   - Motion blur model: σ_blur = ||ω|| * f / (pixel_size * update_rate).
%   - Noise model: σ_i = σ_ref * 10^(0.2*(mag_i - mag_ref)) (SNR scaling).
%   - Weight formula: w ∝ 10^(-0.4*mag) = 1/σ_i^2 (optimal weighting).
%==========================================================================
    fprintf('\n=== Simulating Star Tracker Measurements ===\n');
    
    % Pre-allocate measurement structure array
    meas(N_STR) = createEmptyMeasurement();
    
    % Loop over each star tracker
    for i_str = 1:N_STR
        
        fprintf('\n--- Star Tracker %d ---\n', i_str);
        
        % Compose full transformation: ECI -> Body -> Sensor
        DCM_ECI2Sensor = STR(i_str).DCM_B2S * DCM_true;
        
        % Transform all catalog stars to sensor frame
        r_sensor_all = DCM_ECI2Sensor * catalog.r_ECI;  % [3 x N_stars]
        
        % Filter stars in front of camera (Z > 0) and within FOV cone
        valid_hemisphere = r_sensor_all(3,:) > 0;
        theta_sep_all = acosd(max(min(r_sensor_all(3,:), 1), -1)); % Clamp for stability
        in_fov_logical = valid_hemisphere & (theta_sep_all <= STR(i_str).FOV_deg/2);
        in_fov_idx = find(in_fov_logical);
        
        N_visible = length(in_fov_idx);
        fprintf('  Stars in FOV cone: %d\n', N_visible);
        
        if N_visible < 3
            warning('STR %d: Insufficient stars (%d < 3). Skipping...', i_str, N_visible);
            meas(i_str) = createEmptyMeasurement();
            continue;
        end
        
        % Gnomonic projection to focal plane: (X/Z, Y/Z) * f
        r_sensor = r_sensor_all(:, in_fov_idx);
        x_focal_plane = STR(i_str).focal_length * r_sensor(1,:) ./ r_sensor(3,:); % [m]
        y_focal_plane = STR(i_str).focal_length * r_sensor(2,:) ./ r_sensor(3,:); % [m]
        
        % Convert to pixel coordinates (origin at sensor center)
        pixel_x = x_focal_plane / STR(i_str).pixel_size;
        pixel_y = y_focal_plane / STR(i_str).pixel_size;
        
        % Filter stars within physical sensor boundaries
        half_width  = STR(i_str).resolution(1) / 2;
        half_height = STR(i_str).resolution(2) / 2;
        on_sensor = (abs(pixel_x) <= half_width) & (abs(pixel_y) <= half_height);
        
        pixel_x = pixel_x(on_sensor);
        pixel_y = pixel_y(on_sensor);
        in_fov_idx = in_fov_idx(on_sensor);
        N_visible = length(in_fov_idx);
        
        fprintf('  Stars on sensor:   %d\n', N_visible);
        
        if N_visible < 3
            warning('STR %d: Insufficient stars on sensor (%d < 3). Skipping...', i_str, N_visible);
            meas(i_str) = createEmptyMeasurement();
            continue;
        end
        
        % ===== MODIFICACIÓN: Ruido dependiente de magnitud =====
        % Magnitudes de las estrellas detectadas
        current_mags = catalog.magnitude(in_fov_idx)';
        
        % Magnitud de referencia (estrella más tenue detectable)
        ref_mag = STR(i_str).magnitude_limit;
        
        % Factor de escala de ruido: SNR ∝ sqrt(flux) ∝ 10^(-0.2*mag)
        % Ruido más alto para estrellas tenues, más bajo para brillantes
        noise_scale = 10.^(0.2 * (current_mags - ref_mag));
        
        % Generar ruido Gaussiano con desviación estándar dependiente de magnitud
        noise_x = STR(i_str).centroid_accuracy * noise_scale .* randn(1, N_visible);
        noise_y = STR(i_str).centroid_accuracy * noise_scale .* randn(1, N_visible);
        % ===== FIN MODIFICACIÓN =====
        
        pixel_x_noisy = pixel_x + noise_x;
        pixel_y_noisy = pixel_y + noise_y;
        
        % Add motion blur if angular rate is significant
        motion_blur_std = norm(omega_body) * STR(i_str).focal_length / ...
                          (STR(i_str).pixel_size * STR(i_str).update_rate);
        if motion_blur_std > 0.1
            pixel_x_noisy = pixel_x_noisy + motion_blur_std * randn(1, N_visible);
            pixel_y_noisy = pixel_y_noisy + motion_blur_std * randn(1, N_visible);
            fprintf('  Motion blur:       %.2f pixels (1σ)\n', motion_blur_std);
        end
        
        % Back-project noisy pixel coordinates to unit vectors
        x_focal_noisy = pixel_x_noisy * STR(i_str).pixel_size;
        y_focal_noisy = pixel_y_noisy * STR(i_str).pixel_size;
        z_focal = STR(i_str).focal_length * ones(1, N_visible);
        
        r_sensor_measured = [x_focal_noisy; y_focal_noisy; z_focal];
        r_sensor_measured = r_sensor_measured ./ vecnorm(r_sensor_measured); % Normalize
        
        % Compute flux-based weights (Pogson's law: flux ∝ 10^(-0.4*mag))
        % Ahora consistente con el modelo de ruido: w_i ∝ 1/σ_i^2 ∝ flux
        mag_weights = 10.^(-0.4 * current_mags);
        mag_weights = mag_weights / sum(mag_weights); % Normalize to sum=1
        
        % Store measurements
        meas(i_str).N_stars       = N_visible;
        meas(i_str).star_id       = in_fov_idx;
        meas(i_str).r_body        = STR(i_str).DCM_B2S' * r_sensor_measured; % Sensor -> Body
        meas(i_str).r_ECI_ref     = catalog.r_ECI(:, in_fov_idx);
        meas(i_str).magnitude     = catalog.magnitude(in_fov_idx);
        meas(i_str).weights       = mag_weights';
        meas(i_str).pixel_x_noisy = pixel_x_noisy;
        meas(i_str).pixel_y_noisy = pixel_y_noisy;
        
        fprintf('  Stars identified:  %d\n', N_visible);
    end
end

%==========================================================================
% createEmptyMeasurement: Initialize empty measurement structure.
%==========================================================================
function m = createEmptyMeasurement()
    m.N_stars       = 0;
    m.star_id       = [];
    m.r_body        = [];
    m.r_ECI_ref     = [];
    m.magnitude     = [];
    m.weights       = [];
    m.pixel_x_noisy = [];
    m.pixel_y_noisy = [];
end
