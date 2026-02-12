function IMU = initializeIMU()
%==========================================================================
% initializeIMU - Initialize IMU hardware parameters and error models
%
% OUTPUTS:
%   IMU - Structure with fields:
%         .rate, .dt              - Sampling configuration
%         .DCM_mounting           - Body to IMU mounting misalignment
%         .gyro                   - Gyroscope parameters and error matrices
%         .accel                  - Accelerometer parameters and error matrices
%==========================================================================

    fprintf('\n=== Initializing IMU Model ===\n');
    
    % --- General Configuration ---
    IMU.rate = 120;        % Sampling frequency [Hz]
    IMU.dt   = 1/IMU.rate; % Sampling time step [s]
    
    % Mounting misalignment errors (typical mechanical tolerances)
    rollErr  =  0.002; % Mounting error around X-axis [deg]
    pitchErr = -0.005; % Mounting error around Y-axis [deg]
    yawErr   =  0.003; % Mounting error around Z-axis [deg]
    
    % Construct misalignment DCM (small perturbation from ideal alignment)
    % This matrix is shared between gyro and accel (common mechanical mounting).
    IMU.DCM_mounting = angle2dcm(deg2rad(yawErr), ...
                                 deg2rad(pitchErr), ...
                                 deg2rad(rollErr), ...
                                 'ZYX');
    
    % --- Gyroscope Parameters (Angular Rate) ---
        % 1. Stochastic Noise (Random processes)
            % ARW: High-frequency noise (white noise density). Determines short-term precision.
            IMU.gyro.ARW = 0.05;            % Angle Random Walk [deg/sqrt(h)]
    
            % Bias Instability: Long-term drift (flicker noise 1/f floor).
            IMU.gyro.biasInstability = 1.0; % [deg/h]
    
            % RRW: Rate Random Walk. Models the "wandering" of the bias over time.
            IMU.gyro.RRW = 0.1;             % [deg/h/sqrt(h)] (bias drift rate)
    
        % 2. Deterministic Errors (Calibration residuals)
            % Static Bias: Constant offset per power cycle (Turn-on to turn-on repeatability).
            IMU.gyro.biasStaticLim = 10.0;  % [deg/h] (1σ)
    
            % Scale Factor Error: Sensitivity error (linear gain mismatch).
            IMU.gyro.SF = 300;              % [ppm]
    
            % 3. Geometric Errors
            % Internal Misalignment: Non-orthogonality of the IMU triad axes.
            IMU.gyro.nonortho = 0.05;       % [deg] (1σ)
    
    % --- Accelerometer Parameters (Specific Force) ---
        % 1. Stochastic Noise
        IMU.accel.VRW = 50;                 % Velocity Random Walk [μg/sqrt(Hz)]
        IMU.accel.biasInstability = 20;     % [μg]
        IMU.accel.ARW = 0.05;               % Acceleration Random Walk [m/s/sqrt(h)] (for bias drift)
    
        % 2. Deterministic Errors
        IMU.accel.biasStaticLim = 1.0;      % [mg]
        IMU.accel.SF = 300;                 % [ppm]
     
        % 3. Geometric Errors
        IMU.accel.nonorthogo = 0.05;        % [deg]
    
    % --- Derived Matrices & Pre-Allocation ---
    fprintf('\n=== Initializing IMU Model ===\n');
        
        % --- Gyroscope Deterministic Transform ---
            % Construct Scale Factor Matrix S = diag([sf_x, sf_y, sf_z])
            gyro_SFErr = (IMU.gyro.SF * 1e-6) * randn(3,1); 
            IMU.gyro.M_SF = eye(3) + diag(gyro_SFErr);
            
            fprintf(' Gyro Scale Factor Errors:  [%.1f, %.1f, %.1f] ppm\n', ...
                gyro_SFErr(1)*1e6, gyro_SFErr(2)*1e6, gyro_SFErr(3)*1e6);
            
            % Construct Non-Orthogonality Matrix
            % Upper triangular approximation for small angles.
            gyro_nonorthoErr = deg2rad(IMU.gyro.nonortho) * randn(3,1);
            IMU.gyro.M_nonorth = [                   1, -gyro_nonorthoErr(3),  gyro_nonorthoErr(2);
                                  -gyro_nonorthoErr(3),                    1, -gyro_nonorthoErr(1);
                                   gyro_nonorthoErr(2), -gyro_nonorthoErr(1),                    1];
            
            % Total deterministic transform (Body True -> Gyro IMU Frame)
            IMU.gyro.T_Deterministic = IMU.gyro.M_SF * IMU.gyro.M_nonorth * IMU.DCM_mounting;
            
            % Generate Static Turn-on Bias (Constant for this simulation run)
            IMU.gyro.biasStatic = deg2rad(IMU.gyro.biasStaticLim / 3600) * randn(3,1);
            
            fprintf(' Gyro Static Turn-on Bias:  [%.4f, %.4f, %.4f] deg/h\n', ...
                rad2deg(IMU.gyro.biasStatic(1))*3600, ...
                rad2deg(IMU.gyro.biasStatic(2))*3600, ...
                rad2deg(IMU.gyro.biasStatic(3))*3600);
        
        % --- Accelerometer Deterministic Transform ---
            % Construct Scale Factor Matrix
            accel_SFErr = (IMU.accel.SF * 1e-6) * randn(3,1);
            IMU.accel.M_SF = eye(3) + diag(accel_SFErr);
            
            fprintf(' Accel Scale Factor Errors: [%.1f, %.1f, %.1f] ppm\n', ...
                accel_SFErr(1)*1e6, accel_SFErr(2)*1e6, accel_SFErr(3)*1e6);
            
            % Construct Non-Orthogonality Matrix
            accel_NonorthoErr = deg2rad(IMU.accel.nonorthogo) * randn(3,1);
            IMU.accel.M_nonorth = [                    1, -accel_NonorthoErr(3),  accel_NonorthoErr(2);
                                   -accel_NonorthoErr(3),                     1, -accel_NonorthoErr(1);
                                    accel_NonorthoErr(2), -accel_NonorthoErr(1),                     1];
            
            % Total deterministic transform (Body Frame -> Accel IMU Frame)
            IMU.accel.T_Deterministic = IMU.accel.M_SF * IMU.accel.M_nonorth * IMU.DCM_mounting;
            
            % Generate Static Turn-on Bias
            g0 = 9.80665; % [m/s^2]
            IMU.accel.biasStatic = (IMU.accel.biasStaticLim * 1e-3 * g0) * randn(3,1); % [mg]
            
            fprintf(' Accel Static Turn-on Bias: [%.4f, %.4f, %.4f] mg\n', ...
                IMU.accel.biasStatic(1)/(1e-3*g0), ...
                IMU.accel.biasStatic(2)/(1e-3*g0), ...
                IMU.accel.biasStatic(3)/(1e-3*g0));
    
        fprintf('=== IMU Initialization Complete ===\n');
end
