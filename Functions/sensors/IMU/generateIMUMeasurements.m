function meas = generateIMUMeasurements(t, omegaTrue_body, forceTrue_body, IMU)
%==========================================================================
% generateIMUMeasurements: Generate synthetic IMU measurements (gyro and
% accel) from truth data in body frame using a realistic error model.
%
% Inputs:
%   t                    - Time vector                                         [s], 1xN
%   omegaTrue_body       - True angular rate in body frame                     [rad/s], 3xN
%   forceTrue_body       - True specific force in body frame                   [m/s^2], 3xN
%   IMU                  - IMU parameter structure with required fields:
%                          .gyro.ARW              - Angle Random Walk          [deg/sqrt(h)]
%                          .gyro.RRW              - Rate Random Walk           [deg/h/sqrt(h)]
%                          .gyro.T_Deterministic  - Deterministic transform    3x3
%                          .gyro.biasStatic       - Static (turn-on) gyro bias [rad/s], 3x1
%                          .accel.VRW             - Velocity Random Walk       [µg/sqrt(Hz)]
%                          .accel.ARW             - Bias drift random walk     [m/s/sqrt(h)]
%                          .accel.T_Deterministic - Deterministic transform    3x3
%                          .accel.biasStatic      - Static accel bias          [m/s^2], 3x1
%
% Outputs (struct meas):
%   meas.gyro.omegaBody  - Measured angular rate      [rad/s], 3xN
%   meas.gyro.biasDyn    - Dynamic gyro bias history  [rad/s], 3xN
%   meas.accel.forceBody - Measured specific force    [m/s^2], 3xN
%   meas.accel.biasDyn   - Dynamic accel bias history [m/s^2], 3xN
%==========================================================================

    fprintf('\n=== Simulating IMU Measurements ===\n');

    % Basic checks
    if nargin < 4
        error('generateIMUMeasurements: Not enough input arguments.');
    end

    dt = t(2) - t(1);
    N  = numel(t);

    if size(omegaTrue_body,2) ~= N
        error('omegaTrue_body must have size 3xN with N = length(t).');
    end

    % Pre-allocate outputs
    meas.gyro.omegaBody  = zeros(3, N);
    meas.gyro.biasDyn    = zeros(3, N);
    meas.accel.forceBody = zeros(3, N);
    meas.accel.biasDyn   = zeros(3, N);

    %% ====================================================================
    % 1. GYROSCOPE NOISE & BIAS PARAMETERS (DISCRETE-TIME)
    % =====================================================================

    % 1.1 Angle Random Walk (ARW) -> White noise STD in discrete time
    gyroARW            = deg2rad(IMU.gyro.ARW) / sqrt(3600); % [rad/sqrt(s)]
    gyro_sigmaWhite    = gyroARW / sqrt(dt);                 % [rad/s] per sample

    % 1.2 Rate Random Walk (RRW) -> Bias random walk step
    gyroRRW            = deg2rad(IMU.gyro.RRW) / (3600 * sqrt(3600)); % [rad/s^2/sqrt(s)]
    gyro_sigmaBiasStep = gyroRRW * sqrt(dt);      % [rad/s] per step

    % Initial dynamic bias (start at zero)
    gyroBiasDyn_k = zeros(3,1);

    %% ====================================================================
    % 2. ACCELEROMETER NOISE & BIAS PARAMETERS (DISCRETE-TIME)
    % =====================================================================

    g0 = 9.80665; % [m/s^2]

    % 2.1 Velocity Random Walk (VRW) in [μg/sqrt(Hz)]
    accelVRW         = (IMU.accel.VRW * 1e-6) * g0;   % [m/s^2/sqrt(Hz)]
    accel_sigmaWhite = accelVRW / sqrt(dt);           % [m/s^2] per sample

    % 2.2 Acceleration Random Walk (ARW for bias drift) in [m/s/sqrt(h)]
    accelARW            = IMU.accel.ARW / sqrt(3600); % [m/s^2/sqrt(s)]
    accel_sigmaBiasStep = accelARW * sqrt(dt);        % [m/s^2] per step

    accelBiasDyn_k = zeros(3,1);

    %% ====================================================================
    % 3. MAIN LOOP: APPLY ERROR MODEL SAMPLE BY SAMPLE
    % =====================================================================

    for k = 1:N

        %------------------------------------------------------------------
        % 3.1 TRUE RATE & FORCE IN IMU FRAME (Deterministic transform)
        %------------------------------------------------------------------
        % Apply mounting + scale factor + non-orthogonality:
        omega_IMUTrue =  IMU.gyro.T_Deterministic * omegaTrue_body(:,k);
        force_IMUTrue = IMU.accel.T_Deterministic * forceTrue_body(:,k);

        %------------------------------------------------------------------
        % 3.2 BIAS (Static + Dynamic Random Walk)
        %------------------------------------------------------------------
        % Dynamic bias evolution (RW)
        gyroBiasDyn_k = gyroBiasDyn_k + gyro_sigmaBiasStep * randn(3,1);
        meas.gyro.biasDyn(:,k) = gyroBiasDyn_k;

        accelBiasDyn_k = accelBiasDyn_k + accel_sigmaBiasStep * randn(3,1);
        meas.accel.biasDyn(:,k) = accelBiasDyn_k;

        % Total bias = static (turn-on) + dynamic
         gyroBiasTotal_k =  IMU.gyro.biasStatic + gyroBiasDyn_k;
        accelBiasTotal_k = IMU.accel.biasStatic + accelBiasDyn_k;

        %------------------------------------------------------------------
        % 3.3 WHITE NOISE
        %------------------------------------------------------------------
         gyroNoise_k =  gyro_sigmaWhite * randn(3,1);
        accelNoise_k = accel_sigmaWhite * randn(3,1);

        %------------------------------------------------------------------
        % 3.4 FINAL MEASUREMENT
        %------------------------------------------------------------------
         meas.gyro.omegaBody(:,k) = omega_IMUTrue +  gyroBiasTotal_k +  gyroNoise_k;
        meas.accel.forceBody(:,k) = force_IMUTrue + accelBiasTotal_k + accelNoise_k;

    end
end
