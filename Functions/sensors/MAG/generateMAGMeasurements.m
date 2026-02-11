function mag_meas = generateMAGMeasurements(t, B_true_B, MAG)
%==========================================================================
% generateMAGMeasurements: Generate synthetic magnetometer measurements
%                          from truth magnetic field data using a realistic
%                          error model (Fluxgate/AMR sensor).
%
% Inputs:
%   t                       - Time vector [s], 1xN
%   B_true_B                - True magnetic field vector in body frame [nT], 3xN
%   MAG                     - Magnetometer parameter struct (see mainMAG.m)
%                             Must contain:
%                             .dt, .T_Deterministic, .biasHardIron, .sigma_white
%                             .biasInstability, .range, .resolution
%
% Outputs (struct mag_meas):
%   mag_meas.B_meas         - Measured magnetic field [nT], 3xN
%   mag_meas.bias_hard_dyn  - Dynamic bias history (drift) [nT], 3xN
%   mag_meas.B_clean        - Signal before quantization/noise (debug)
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
%==========================================================================

    % --- Basic Input Validation ---
    if nargin < 3
        error('generateMAGMeasurements: Not enough input arguments.');
    end

    N = numel(t);
    if ~isequal(size(B_true_B), [3, N])
        error('B_true_B must have size 3xN with N = length(t).');
    end

    % Pre-allocate outputs
    mag_meas.B_meas     = zeros(3, N);
    mag_meas.B_clean    = zeros(3, N);
    mag_meas.B_true_S   = zeros(3, N);
    mag_meas.B_det      = zeros(3, N);
    mag_meas.bias_dyn   = zeros(3, N);
    mag_meas.bias_total = zeros(3, N);

    % Noise (white)
    sigmaWhite = MAG.sigmaWhite;

    % Bounded bias drift (1st-order Gauss-Markov)
    tau = MAG.biasTau;
    if tau <= 0
        a = 0;
    else
        a = exp(-MAG.dt / tau);
    end
    sigmaSS = MAG.biasInstability;       % Steady-State 1-sigma (per axis)
    q = sqrt(max(1 - a^2, 0)) * sigmaSS; % Driving noise std per step

    bias_dyn = zeros(3,1);

    % Limits
    max_range  = MAG.range;
    resolution = MAG.resolution;

    for k = 1:N
        % 1) Body truth -> Sensor truth (mounting applied)
        B_true_S = MAG.DCM_B2S_true * B_true_B(:,k);
        mag_meas.B_true_S(:,k) = B_true_S;

        % 2) Deterministic distortion in sensor axes
        B_det = MAG.M_Deterministic * B_true_S;
        mag_meas.B_det(:,k) = B_det;

        % 3) Dynamic bias update (Gauss-Markov)
        bias_dyn = a * bias_dyn + q * randn(3,1);
        mag_meas.bias_dyn(:,k) = bias_dyn;

        % 4) Total bias (static + dynamic)
        b_total = MAG.biasHardIron + bias_dyn;
        mag_meas.bias_total(:,k) = b_total;

        % 5) White noise
        noise_k = sigmaWhite * randn(3,1);

        % 6) Analog (pre-quantization)
        B_analog = B_det + b_total + noise_k;
        mag_meas.B_clean(:,k) = B_analog;

        % 7) Saturation
        B_clip = max(min(B_analog, max_range), -max_range);

        % 8) Quantization
        if resolution > 0
            B_dig = round(B_clip / resolution) * resolution;
        else
            B_dig = B_clip;
        end

        mag_meas.B_meas(:,k) = B_dig;
    end
end
