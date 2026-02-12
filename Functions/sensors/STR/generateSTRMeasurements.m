function meas = generateSTRMeasurements(STR, nSTR, catalog, DCMTrue, omegaBody)
%==========================================================================
% generateSTRMeasurements: Simulate star tracker observations with
%                          realistic noise models.
%
% Inputs:
%   STR         - Structure array (1 x nSTR) with STR specifications:
%                 .DCM_B2S           - Body-to-Sensor DCM (3x3)
%                 .FOV               - Field of view [deg]
%                 .resolution        - Sensor dimensions [pixels]
%                 .pixelSize         - Pixel pitch [m]
%                 .focalLength       - Optical focal length [m]
%                 .centroidAccuracy  - Centroiding std dev [pixels] for 
%                                      reference magnitude
%                 .magnitudeLimit   - Reference magnitude for centroidAccuracy
%                 .rate       - Measurement frequency [Hz]
%   nSTR        - Number of active star trackers (1 or 2).
%   catalog     - Star catalog structure:
%                 .rECI             - 3xN unit vectors in ECI frame
%                 .magnitude         - Nx1 absolute magnitudes
%   DCM_true    - True ECI-to-Body DCM (3x3).
%   omega_body  - Satellite angular velocity in body frame [rad/s].
%
% Outputs:
%   meas        - Structure array (1 x nSTR) with measurements:
%                 .nStars           - Number of detected stars
%                 .starID           - Catalog indices of detected stars
%                 .rBody            - 3xM unit vectors in body frame
%                 .rECI_ref         - 3xM reference unit vectors (ECI)
%                 .magnitude         - Mx1 absolute magnitudes
%                 .weights           - Mx1 measurement weights (flux-based)
%                 .pixelX_noisy     - Mx1 pixel coordinates (for viz)
%                 .pixelY_noisy     - Mx1 pixel coordinates (for viz)
%
% Method:
%   1. Transform catalog stars to each sensor frame via DCM_B2S * DCMTrue.
%   2. Apply gnomonic projection to focal plane (pinhole camera model).
%   3. Filter stars within FOV and physical sensor boundaries.
%   4. Add magnitude-dependent Gaussian centroiding noise and motion blur.
%   5. Back-project noisy pixel coordinates to unit vectors.
%   6. Compute flux-based weights (brighter stars = higher weight).
%
% Notes:
%   - Minimum 3 stars required per STR for attitude determination.
%   - Motion blur model: σ_blur = ||ω|| * f / (pixelSize * rate).
%   - Noise model: σ_i = σ_ref * 10^(0.2*(mag_i - mag_ref)) (SNR scaling).
%   - Weight formula: w ∝ 10^(-0.4*mag) = 1/σ_i^2 (optimal weighting).
%==========================================================================

    fprintf('\n=== Simulating Star Tracker Measurements ===\n');
    
    % Pre-allocate measurement structure array
    meas(nSTR) = createEmptySTRMeasurement();
    
    % Loop over each star tracker
    for iSTR = 1:nSTR
        
        fprintf('\n--- Star Tracker %d ---\n', iSTR);
        
        % Compose full transformation: ECI -> Body -> Sensor
        DCM_ECI2STR = STR(iSTR).DCM_B2S * DCMTrue;
        
        % Transform all catalog stars to sensor frame
        rSTR_all = DCM_ECI2STR * catalog.rECI;  % [3 x nStars]
        
        % Filter stars in front of camera (Z > 0) and within FOV cone
        validHemisphere   = rSTR_all(3,:) > 0;
        separationAng_all = acosd(max(min(rSTR_all(3, :), 1), -1)); % Clamp for stability
        inFOV_logical     = validHemisphere & (separationAng_all <= STR(iSTR).FOV/2);
        inFOV_idx         = find(inFOV_logical);
        
        nVisible = length(inFOV_idx);
        fprintf('  Stars in FOV cone: %d\n', nVisible);
        
        if nVisible < 3
            warning('STR %d: Insufficient stars (%d < 3). Skipping...', iSTR, nVisible);
            meas(iSTR) = createEmptySTRMeasurement();
            continue;
        end
        
        % Gnomonic projection to focal plane: (X/Z, Y/Z) * f
        rSTR = rSTR_all(:, inFOV_idx);
        x_focalPlane = STR(iSTR).focalLength * rSTR(1,:) ./ rSTR(3,:); % [m]
        y_focalPlane = STR(iSTR).focalLength * rSTR(2,:) ./ rSTR(3,:); % [m]
        
        % Convert to pixel coordinates (origin at sensor center)
        pixelX = x_focalPlane / STR(iSTR).pixelSize;
        pixelY = y_focalPlane / STR(iSTR).pixelSize;
        
        % Filter stars within physical sensor boundaries
        half_width  = STR(iSTR).resolution(1) / 2;
        half_height = STR(iSTR).resolution(2) / 2;
        on_sensor = (abs(pixelX) <= half_width) & (abs(pixelY) <= half_height);
        
        pixelX = pixelX(on_sensor);
        pixelY = pixelY(on_sensor);
        inFOV_idx = inFOV_idx(on_sensor);
        nVisible = length(inFOV_idx);
        
        fprintf('  Stars on sensor:   %d\n', nVisible);
        
        if nVisible < 3
            warning('STR %d: Insufficient stars on sensor (%d < 3). Skipping...', iSTR, nVisible);
            meas(iSTR) = createEmptySTRMeasurement();
            continue;
        end
        
        % Magnitudes of the detected stars
        currentMags = catalog.magnitude(inFOV_idx)';
        
        % Reference magnitude (faintest detectable star)
        refMag = STR(iSTR).magnitudeLimit;
        
        % Noise scaling factor: SNR ∝ sqrt(flux) ∝ 10^(-0.2*mag)
        % Higher noise for fainter stars, lower for brighter ones
        noiseScale = 10.^(0.2 * (currentMags - refMag));
        
        % Generate Gaussian noise with magnitude-dependent standard deviation
        noiseX = STR(iSTR).centroidAccuracy * noiseScale .* randn(1, nVisible);
        noiseY = STR(iSTR).centroidAccuracy * noiseScale .* randn(1, nVisible);
        
        pixelX_noisy = pixelX + noiseX;
        pixelY_noisy = pixelY + noiseY;
        
        % Add motion blur if angular rate is significant
        motionBlur_std = norm(omegaBody) * STR(iSTR).focalLength / ...
                         (STR(iSTR).pixelSize * STR(iSTR).rate);
        if motionBlur_std > 0.1
            pixelX_noisy = pixelX_noisy + motionBlur_std * randn(1, nVisible);
            pixelY_noisy = pixelY_noisy + motionBlur_std * randn(1, nVisible);
            fprintf('  Motion blur:       %.2f pixels (1σ)\n', motionBlur_std);
        end
        
        % Back-project noisy pixel coordinates to unit vectors
        x_focalNoisy = pixelX_noisy * STR(iSTR).pixelSize;
        y_focalNoisy = pixelY_noisy * STR(iSTR).pixelSize;
        z_focal = STR(iSTR).focalLength * ones(1, nVisible);
        
        r_STRMeasured = [x_focalNoisy; y_focalNoisy; z_focal];
        r_STRMeasured = r_STRMeasured ./ vecnorm(r_STRMeasured); % Normalize
        
        % Compute flux-based weights (Pogson's law: flux ∝ 10^(-0.4*mag))
        magWeights = 10.^(-0.4 * currentMags);
        magWeights = magWeights / sum(magWeights); % Normalize to sum=1
        
        % Store measurements
        meas(iSTR).nStars       = nVisible;
        meas(iSTR).starID       = inFOV_idx;
        meas(iSTR).rBody        = STR(iSTR).DCM_B2S' * r_STRMeasured; % Sensor -> Body
        meas(iSTR).rECI_ref     = catalog.rECI(:, inFOV_idx);
        meas(iSTR).magnitude    = catalog.magnitude(inFOV_idx);
        meas(iSTR).weights      = magWeights';
        meas(iSTR).pixelX_noisy = pixelX_noisy;
        meas(iSTR).pixelY_noisy = pixelY_noisy;
        
        fprintf('  Stars identified:  %d\n', nVisible);
    end
end

%==========================================================================
% createEmptyMeasurement: Initialize empty STR measurement structure.
%==========================================================================
function m = createEmptySTRMeasurement()
    m.nStars       = 0;
    m.starID       = [];
    m.rBody        = [];
    m.rECI_ref     = [];
    m.magnitude    = [];
    m.weights      = [];
    m.pixelX_noisy = [];
    m.pixelY_noisy = [];
end
