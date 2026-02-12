function ekf = initializeMEKF(IMU, MAG, STR)
%==========================================================================
% initializeMEKF - Initialize Multiplicative Extended Kalman Filter (MEKF)
%                  for spacecraft attitude estimation
%
% INPUTS:
%   IMU - IMU parameter structure (from initializeIMU.m)
%   MAG - Magnetometer parameter structure (from initializeMAG.m)
%   STR - Star Tracker parameter array (from initializeSTR.m)
%
% OUTPUTS:
%   ekf - EKF state structure with fields:
%         .x      - Error state vector [6x1]: [δθ; b_gyro]
%         .P      - Error covariance matrix [6x6]
%         .q_nom  - Nominal attitude quaternion [4x1] (scalar-first)
%         .Q      - Process noise covariance [6x6]
%         .R_mag  - Magnetometer measurement noise [3x3]
%         .R_str  - Star Tracker measurement noise [3x3]
%
% MEKF STATE DEFINITION:
%   x = [δθ (3x1); b_gyro (3x1)]
%   where:
%     δθ     : Attitude error (small angle, Gibbs vector or MRP)
%     b_gyro : Gyroscope bias error [rad/s]
%
%   The true attitude is maintained separately as a quaternion:
%     q_true = δq(δθ) ⊗ q_nom
%==========================================================================

    fprintf('\n=== Initializing MEKF ===\n');
    
    %% ===================================================================
    % 1. ERROR STATE INITIALIZATION
    % ====================================================================
    % Initial error state (zero error assumption)
    ekf.x = zeros(6,1);  % [δθ; b_gyro]
    
    fprintf(' Error state dimension: %d\n', length(ekf.x));
    fprintf('   - Attitude error (δθ):   3 states\n');
    fprintf('   - Gyro bias error:       3 states\n');
    
    %% ===================================================================
    % 2. INITIAL COVARIANCE MATRIX
    % ====================================================================
    % Attitude uncertainty: assume ±10° initial pointing error (1σ)
    sigma_att0 = deg2rad(10);  % [rad]
    P_att      = (sigma_att0^2) * eye(3);
    
    % Gyro bias uncertainty: use turn-on bias repeatability spec
    sigma_bias0 = deg2rad(IMU.gyro.biasStaticLim / 3600);  % [deg/h] -> [rad/s]
    P_bias      = (sigma_bias0^2) * eye(3);
    
    ekf.P = blkdiag(P_att, P_bias);
    
    fprintf(' Initial covariance:\n');
    fprintf('   - Attitude uncertainty:  %.2f deg (1σ)\n', rad2deg(sqrt(P_att(1,1))));
    fprintf('   - Gyro bias uncertainty: %.4f deg/h (1σ)\n', ...
            rad2deg(sqrt(P_bias(1,1))) * 3600);
    
    %% ===================================================================
    % 3. NOMINAL ATTITUDE (will be updated with TRIAD or initial guess)
    % ====================================================================
    ekf.q_nom = [1; 0; 0; 0];  % Identity quaternion (ECI = Body initially)
    
    fprintf(' Nominal attitude:        Identity (will be initialized with TRIAD)\n');
    
    %% ===================================================================
    % 4. PROCESS NOISE COVARIANCE (from Allan Deviation analysis)
    % ====================================================================
    ekf.Q = computeProcessNoise(IMU);
    
    fprintf(' Process noise (Q):\n');
    fprintf('   - ARW noise (gyro):      %.2e (rad/s)²/s\n', ekf.Q(1,1));
    fprintf('   - RRW noise (bias drift):%.2e (rad/s)²/s\n', ekf.Q(4,4));
    
    %% ===================================================================
    % 5. MEASUREMENT NOISE COVARIANCES
    % ====================================================================
    % Magnetometer noise
    ekf.R_mag = computeMagNoise(MAG);
    fprintf(' Magnetometer noise (R):  %.2f nT (1σ per axis)\n', sqrt(ekf.R_mag(1,1)));
    
    % Star Tracker noise
    ekf.R_str = computeSTRNoise(STR);
    fprintf(' Star Tracker noise (R):  %.2f arcsec (1σ cross-boresight)\n', ...
            rad2deg(sqrt(ekf.R_str(1,1))) * 3600);
    
    fprintf('=== MEKF Initialization Complete ===\n');
end

%% =======================================================================
% HELPER FUNCTION: Compute Process Noise from IMU Allan Deviation Parameters
% ========================================================================
function Q = computeProcessNoise(IMU)
    % Process noise for MEKF with state [δθ; b_gyro]
    %
    % The continuous-time process noise spectral density is:
    %   Q_c = blkdiag(σ_ARW² * I₃, σ_RRW² * I₃)
    %
    % Discrete-time approximation (for small dt):
    %   Q_d ≈ Q_c * dt
    
    %% 1. Attitude Error Process Noise (from Angle Random Walk)
    % ARW represents the high-frequency white noise in gyro measurements
    % Units: [deg/sqrt(h)] -> [rad/sqrt(s)]
    ARW_rad_sqrtS = deg2rad(IMU.gyro.ARW) / sqrt(3600);  % [rad/√s]
    
    % Spectral density for attitude propagation error
    % δθ is driven by gyro white noise: σ_δθ = ARW * √dt
    sigma_ARW = ARW_rad_sqrtS^2;  % [(rad/s)²/Hz] = [(rad/s)²·s]
    
    Q_att = sigma_ARW * eye(3);
    
    %% 2. Gyro Bias Process Noise (from Rate Random Walk)
    % RRW represents the drift rate of the gyro bias
    % Units: [deg/h/sqrt(h)] -> [rad/s²]
    RRW_rad_s2 = deg2rad(IMU.gyro.RRW) / 3600 / sqrt(3600);  % [rad/s²]
    
    % Spectral density for bias random walk
    sigma_RRW = RRW_rad_s2^2;  % [(rad/s²)²/Hz] = [(rad/s)²/s]
    
    Q_bias = sigma_RRW * eye(3);
    
    %% 3. Assemble Block-Diagonal Process Noise Matrix
    Q = blkdiag(Q_att, Q_bias);
end

%% =======================================================================
% HELPER FUNCTION: Compute Magnetometer Measurement Noise
% ========================================================================
function R_mag = computeMagNoise(MAG)
    % Magnetometer measurement noise covariance
    %
    % The magnetometer provides a 3-axis vector measurement of the 
    % magnetic field. The noise is assumed white, Gaussian, and isotropic.
    %
    % Noise sources:
    %   1. White noise (σ_white from noise density * √bandwidth)
    %   2. Bias instability (modeled as process noise, not measurement)
    %   3. Quantization noise (negligible for modern ADCs)
    
    % White noise RMS (already computed in initializeMAG)
    sigma_white_nT = MAG.sigmaWhite;  % [nT]
    
    % Measurement noise covariance (isotropic, uncorrelated axes)
    R_mag = (sigma_white_nT^2) * eye(3);  % [nT²]
    
    % Note: We use the raw field measurement, not normalized direction.
    % For direction-only measurements (unit vectors), R would be scaled
    % by 1/||B||² to account for normalization.
end

%% =======================================================================
% HELPER FUNCTION: Compute Star Tracker Measurement Noise
% ========================================================================
function R_str = computeSTRNoise(STR)
    % Star Tracker measurement noise covariance
    %
    % STR provides attitude quaternions or star unit vectors. The noise
    % is determined by:
    %   1. Centroiding error (dominant for short exposures)
    %   2. Motion blur (becomes significant above ~0.1 deg/s)
    %   3. Catalog errors (typically negligible)
    
    %% 1. Single-Star Centroiding Noise
    % Convert pixel error to angular error
    sigma_centroid_pix = STR(1).centroidAccuracy;  % [pixels, 1σ]
    
    % Angular noise per star [rad]
    sigma_star_rad = sigma_centroid_pix * STR(1).pixelSize / STR(1).focalLength;
    
    %% 2. Attitude Noise (assumes ~10 stars per frame, typical)
    % Multi-star averaging improves accuracy: σ_att ≈ σ_star / √N
    N_stars_avg = 10;  % Conservative estimate
    
    % Cross-boresight accuracy (better constrained)
    sigma_crossbore = sigma_star_rad / sqrt(N_stars_avg);
    
    % Roll-axis accuracy (worse due to cylindrical symmetry)
    % Typically 3-5x worse than cross-boresight
    sigma_roll = sigma_crossbore * 4;
    
    %% 3. Measurement Noise Covariance (Attitude Error, Body Frame)
    % For MEKF, we model the innovation as attitude error δθ
    % The measurement is effectively a 3-axis attitude correction
    
    % Option A: Isotropic (conservative, simple)
    R_str = (sigma_crossbore^2) * eye(3);
    
    % Option B: Anisotropic (more accurate, requires knowledge of STR boresight)
    % Assume STR boresight is +Z body axis
    % R_str = diag([sigma_crossbore^2, sigma_crossbore^2, sigma_roll^2]);
    
    % We use Option A for simplicity. Anisotropy can be added if STR
    % mounting is known and performance is critical.
end
