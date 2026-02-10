function imu_meas = generateIMUMeasurements(t, omega_true_b, f_true_b, IMU)
%==========================================================================
% generateIMUMeasurements: Generate synthetic IMU measurements (gyro and
% accelerometer) from truth data in body frame using a realistic error model.
%
% Inputs:
%   t                       - Time vector [s], 1xN
%   omega_true_b            - True angular rate in BODY frame [rad/s], 3xN
%   f_true_b                - True specific force in BODY frame [m/s^2], 3xN
%                             (set [] if accelerometer is not used)
%   IMU                     - IMU parameter struct (see mainIMU.m)
%
% Outputs (struct imu_meas):
%   imu_meas.gyro_meas_b    - Measured angular rate [rad/s], 3xN
%   imu_meas.gyro_bias_dyn  - Dynamic gyro bias history [rad/s], 3xN
%   imu_meas.accel_meas_b   - Measured specific force [m/s^2], 3xN (if used)
%   imu_meas.accel_bias_dyn - Dynamic accel bias history [m/s^2], 3xN
%
% Notes:
%   - Truth inputs (omega_true_b, f_true_b) are assumed already expressed
%     in the BODY frame. Orbit, gravity model, etc., are handled upstream
%     by a dynamics/trajectory generator, NOT inside this function.
%   - Error model includes:
%       * Static bias (turn-on)
%       * Dynamic bias (random walk from Bias Instability / RRW)
%       * Scale factor errors
%       * Axes misalignment / non-orthogonality
%       * White noise (ARW / VRW)
%==========================================================================

    % Basic checks
    if nargin < 4
        error('generateIMUMeasurements: Not enough input arguments.');
    end

    dt = t(2) - t(1);
    N  = numel(t);

    if size(omega_true_b,2) ~= N
        error('omega_true_b must have size 3xN with N = length(t).');
    end

    % Pre-allocate outputs
    imu_meas.gyro_meas_b   = zeros(3, N);
    imu_meas.gyro_bias_dyn = zeros(3, N);
    imu_meas.accel_meas_b   = zeros(3, N);
    imu_meas.accel_bias_dyn = zeros(3, N);

    %% ====================================================================
    % 1. GYROSCOPE NOISE & BIAS PARAMETERS (DISCRETE-TIME)
    % =====================================================================

    % 1.1 Angle Random Walk (ARW) -> White noise std in discrete time
    % ARW is given in [deg/sqrt(h)] in IMU.gyro.ARW.
    ARW_deg_sqrt_hr    = IMU.gyro.ARW;
    ARW_rad_sqrt_s     = deg2rad(ARW_deg_sqrt_hr) / sqrt(3600);             % [rad/sqrt(s)]
    sigma_gyro_white   = ARW_rad_sqrt_s / sqrt(dt);                         % [rad/s] per sample

    % 1.2 Rate Random Walk (RRW) -> Bias random walk step
    % RRW given approx in [deg/h/sqrt(h)].
    RRW_deg_hr_sqrt_hr = IMU.gyro.RRW;
    RRW_rad_s_sqrt_s   = deg2rad(RRW_deg_hr_sqrt_hr) / (3600 * sqrt(3600)); % [rad/s^2/sqrt(s)]
    sigma_bias_step    = RRW_rad_s_sqrt_s * sqrt(dt);                       % [rad/s] per step

    % Initial dynamic bias (start at zero; static bias is IMU.gyro.biasStatic)
    bias_dyn_k = zeros(3,1);

    %% ====================================================================
    % 2. ACCELEROMETER NOISE & BIAS PARAMETERS (DISCRETE-TIME)
    % =====================================================================
    g0 = 9.80665; % [m/s^2]

    % 2.1 Velocity Random Walk (VRW) in [micro-g/sqrt(Hz)]
    VRW_ug_sqrt_Hz    = IMU.accel.VRW;
    VRW_mps2_sqrt_Hz  = (VRW_ug_sqrt_Hz * 1e-6) * g0;       % [m/s^2/sqrt(Hz)]
    sigma_accel_white = VRW_mps2_sqrt_Hz / sqrt(dt);        % [m/s^2] per sample

    % 2.2 Acceleration Random Walk (ARW for bias drift) in [m/s/sqrt(h)]
    ARW_accel_ms_sqrt_hr = IMU.accel.ARW;
    ARW_accel_ms2_sqrt_s = ARW_accel_ms_sqrt_hr / sqrt(3600); % [m/s^2/sqrt(s)]
    sigma_accel_bias_step = ARW_accel_ms2_sqrt_s * sqrt(dt);  % [m/s^2] per step

    bias_accel_dyn_k = zeros(3,1);

    %% ====================================================================
    % 3. MAIN LOOP: APPLY ERROR MODEL SAMPLE BY SAMPLE
    % =====================================================================

    for k = 1:N

        %------------------------------------------------------------------
        % 3.1 TRUE RATES IN SENSOR FRAME (Deterministic transform)
        %------------------------------------------------------------------
        % Apply mounting + scale factor + non-orthogonality:
        % omega_sensor_true = T_deterministic * omega_true_b
        omega_sensor_true = IMU.gyro.T_Deterministic * omega_true_b(:,k);

        %------------------------------------------------------------------
        % 3.2 GYRO BIAS (STATIC + DYNAMIC RANDOM WALK)
        %------------------------------------------------------------------
        % Dynamic bias evolution (random walk)
        bias_dyn_k = bias_dyn_k + sigma_bias_step * randn(3,1);
        imu_meas.gyro_bias_dyn(:,k) = bias_dyn_k;

        % Total bias = static (turn-on) + dynamic
        bias_total_k = IMU.gyro.biasStatic + bias_dyn_k;

        %------------------------------------------------------------------
        % 3.3 GYRO WHITE NOISE
        %------------------------------------------------------------------
        noise_gyro_k = sigma_gyro_white * randn(3,1);

        %------------------------------------------------------------------
        % 3.4 FINAL GYRO MEASUREMENT
        %------------------------------------------------------------------
        imu_meas.gyro_meas_b(:,k) = omega_sensor_true + bias_total_k + noise_gyro_k;

        %------------------------------------------------------------------
        % 3.5 ACCELEROMETER MODEL
        %------------------------------------------------------------------
        % Deterministic transform for accelerometer triad
        f_sensor_true = IMU.accel.T_Deterministic * f_true_b(:,k);

        % Dynamic bias random walk
        bias_accel_dyn_k = bias_accel_dyn_k + sigma_accel_bias_step * randn(3,1);
        imu_meas.accel_bias_dyn(:,k) = bias_accel_dyn_k;

        % Static accel bias (turn-on), modeled similar to gyro
        if ~isfield(IMU.accel, 'biasStatic')
            % Convert [mg] to [m/s^2]
            b_lim = IMU.accel.biasStaticLimit * 1e-3 * g0;  % [mg] -> [m/s^2]
            IMU.accel.biasStatic = b_lim * randn(3,1);
        end

        bias_accel_total_k = IMU.accel.biasStatic + bias_accel_dyn_k;

        % White noise
        noise_accel_k = sigma_accel_white * randn(3,1);

        % Final accel measurement
        imu_meas.accel_meas_b(:,k) = f_sensor_true + bias_accel_total_k + noise_accel_k;
    end
end
