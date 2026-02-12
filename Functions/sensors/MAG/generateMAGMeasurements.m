function meas = generateMAGMeasurements(t, BTrue_body, MAG)
%==========================================================================
% generateMAGMeasurements: Generate synthetic magnetometer measurements
%                          from truth magnetic field data using a realistic
%                          error model (Fluxgate/AMR MAG).
%
% Inputs:
%   t                  - Time vector                                    [s], 1xN
%   BTrue_body         - True magnetic field vector in body frame       [nT], 3xN
%   MAG                - Magnetometer parameter structure with required fields:
%                        .dt              - Sampling time               [s]
%                        .M_Deterministic - Deterministic transform     3x3
%                        .biasHardIron    - Static hard iron bias       [nT], 3x1
%                        .sigmaWhite      - White noise RMS             [nT]
%                        .biasTau         - Bias correlation time       [s]
%                        .biasInstability - Steady-state bias sigma     [nT]
%                        .range           - Full scale range            [nT]
%                        .resolution      - Digital quantization step   [nT]
%
% Outputs (struct meas):
%   meas.B_body        - Measured magnetic field (Digital)              [nT], 3xN
%   meas.biasDyn       - Dynamic bias history (Drift)                   [nT], 3xN
%   meas.BClean        - Analog signal (Deterministic + Bias + Noise)   [nT], 3xN
%
% Algorithm:
%   The measurement model simulates the signal path of a vector magnetometer:
%     B_meas = Quantize( Saturate( T_det * B_true + b_hard + b_dyn + noise ) )
%
%   1. Deterministic Distortion (Soft Iron, Non-Orthogonality, Scale Factor):
%      B_det = M_Deterministic * B_true_body
%
%   2. Bias Application:
%      - Static Hard Iron: Constant offset (biasHardIron)
%      - Dynamic Drift: 1st-order Gauss-Markov process (bDyn)
%        bDyn[k+1] = exp(-dt/tau)*bDyn[k] + wBias[k]
%
%   3. MAG Noise:
%      - Add wideband white Gaussian noise (sigmaWhite)
%
%   4. Analog Signal Conditioning:
%      B_analog = BDet + biasHardIron + bDyn + noise
%
%   5. Digitization (ADC):
%      - Saturation to full-scale range (+/- MAG.range)
%      - Uniform quantization (MAG.resolution)
%==========================================================================

    fprintf('\n=== Simulating IMU Measurements ===\n');

    % -Basic input validation
    if nargin < 3
        error('generateMAGMeasurements: Not enough input arguments.');
    end

    N = numel(t);
    if ~isequal(size(BTrue_body), [3, N])
        error('BTrue_body must have size 3xN with N = length(t).');
    end

    % Pre-allocate outputs
    meas.B         = zeros(3, N);
    meas.BClean    = zeros(3, N);
    meas.BDet      = zeros(3, N);
    meas.biasDyn   = zeros(3, N);
    meas.biasTotal = zeros(3, N);

    % Whie noise
    sigmaWhite = MAG.sigmaWhite;

    % Bounded bias drift (1st-order Gauss-Markov)
    tau = MAG.biasTau;
    if tau <= 0
        a = 0;
    else
        a = exp(-MAG.dt / tau);
    end
    sigmaSS = MAG.biasInstability;             % Steady-State 1σ (per axis)
    q       = sqrt(max(1 - a^2, 0)) * sigmaSS; % Driving noise STD per step

    biasDyn = zeros(3,1);

    % Limits
    maxRange   = MAG.range;
    resolution = MAG.resolution;

    for k = 1:N
        % 1) Deterministic distortion in MAG axes
        BDet = MAG.M_Deterministic * BTrue_body(:,k);
        meas.BDet(:,k) = BDet;

        % 2) Dynamic bias update (Gauss-Markov)
        biasDyn = a * biasDyn + q * randn(3,1);
        meas.biasDyn(:,k) = biasDyn;

        % 3) Total bias (static + dynamic)
        biasTotal = MAG.biasHardIron + biasDyn;
        meas.biasTotal(:,k) = biasTotal;

        % 4) White noise
        noise_k = sigmaWhite * randn(3,1);

        % 5) Analog (pre-quantization)
        BAnalog = BDet + biasTotal + noise_k;
        meas.BClean(:,k) = BAnalog;

        % 6) Saturation
        BClip = max(min(BAnalog, maxRange), -maxRange);

        % 7) Quantization
        if resolution > 0
            BDig = round(BClip / resolution) * resolution;
        else
            BDig = BClip;
        end

        meas.B(:,k) = BDig;
    end
end
