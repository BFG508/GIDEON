function mag_meas = generateMAGMeasurements(t, B_true_B_nT, MAG)
%==========================================================================
% generateMAGMeasurements: Generate synthetic Magnetometer measurements
%                          from truth magnetic field data using a realistic
%                          error model (Fluxgate/AMR sensor).
%
% Inputs:
%   t            - Time vector [s], 1xN
%   B_true_B_nT  - True magnetic field vector in BODY frame [nT], 3xN
%   MAG          - Magnetometer parameter struct (see mainMAG.m)
%                  Must contain:
%                  .dt, .T_Deterministic, .biasHardIron, .sigma_white
%                  .biasInstability, .range, .resolution
%
% Outputs (struct mag_meas):
%   mag_meas.B_meas_nT      - Measured magnetic field [nT], 3xN
%   mag_meas.bias_hard_dyn  - Dynamic bias history (drift) [nT], 3xN
%   mag_meas.B_clean_nT     - Signal before quantization/noise (debug)
%
% Algorithm:
%   B_meas = Quantize( T_det * B_true + b_hard_static + b_dynamic + noise )
%
%   1. Deterministic Transform (Soft Iron + Non-Orthogonality + Mounting)
%   2. Add Static Hard Iron Bias
%   3. Add Dynamic Bias (Random Walk / Flicker model)
%   4. Add Wideband Gaussian Noise
%   5. Saturate to Dynamic Range
%   6. Quantize to Digital Resolution
%
%==========================================================================

    % --- Basic Input Validation ---
    if nargin < 3
        error('generateMAGMeasurements: Not enough input arguments.');
    end
    
    N = numel(t);
    if size(B_true_B_nT, 2) ~= N
        error('B_true_B_nT must have size 3xN with N = length(t).');
    end

    % --- Pre-allocate Outputs ---
    mag_meas.B_meas_nT      = zeros(3, N);
    mag_meas.bias_hard_dyn  = zeros(3, N);
    mag_meas.B_clean_nT     = zeros(3, N); % Useful for debugging
    
    % --- Noise Model Parameters ---
    % 1. Wideband Noise (White)
    % MAG.sigma_white is typically already scaled for the sample rate in mainMAG
    sigma_white = MAG.sigma_white; 
    
    % 2. Bias Instability (Modeled as Random Walk)
    % Bias Instability is usually given as 1-sigma variation over a period.
    % We model it as a random walk step to simulate "drift".
    % Approximation: step_std = sigma_instability * sqrt(dt / T_correlation)
    % For generic sim, we use a small step derived from instability.
    % Heuristic: drift reaches 1-sigma in approx 1 hour (3600s)
    bias_instability_nT = MAG.biasInstability; % [nT]
    sigma_bias_step = bias_instability_nT * sqrt(MAG.dt / 3600); % [nT/step]
    
    % Initialize dynamic bias state
    bias_dyn_k = zeros(3,1);
    
    % --- Hardware Limits ---
    max_range  = MAG.range;      % [nT]
    resolution = MAG.resolution; % [nT]

    fprintf('   Simulating %d magnetometer samples...\n', N);

    %% ====================================================================
    %  MAIN SIMULATION LOOP
    % =====================================================================
    for k = 1:N
        
        % 1. Get True Field Vector (Body Frame)
        B_true_k = B_true_B_nT(:, k);
        
        % 2. Apply Deterministic Errors
        % B_det = (M_Soft * M_NonOrth * DCM_Mount) * B_true
        B_det = MAG.T_Deterministic * B_true_k;
        
        % 3. Update Dynamic Bias (Random Walk)
        bias_dyn_k = bias_dyn_k + sigma_bias_step * randn(3,1);
        mag_meas.bias_hard_dyn(:, k) = bias_dyn_k;
        
        % 4. Total Bias (Static Hard Iron + Dynamic Drift)
        b_total = MAG.biasHardIron + bias_dyn_k;
        
        % 5. Add Stochastic Noise (White Gaussian)
        noise_k = sigma_white * randn(3,1);
        
        % 6. Compile Analog Signal (Infinite Precision)
        B_analog = B_det + b_total + noise_k;
        mag_meas.B_clean_nT(:, k) = B_analog; % Store pre-quantization
        
        % 7. Sensor Saturation (Clipping)
        B_clipped = max(min(B_analog, max_range), -max_range);
        
        % 8. ADC Quantization (Digital Resolution)
        if resolution > 0
            B_digital = round(B_clipped / resolution) * resolution;
        else
            B_digital = B_clipped;
        end
        
        % Store final measurement
        mag_meas.B_meas_nT(:, k) = B_digital;
        
    end
    
    fprintf('   Done.\n');

end
