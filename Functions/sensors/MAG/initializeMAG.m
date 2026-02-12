function MAG = initializeMAG()
%==========================================================================
% initializeMAG - Initialize magnetometer hardware parameters and error models
%
% OUTPUTS:
%   MAG - Structure with fields:
%         .rate, .dt              - Sampling configuration
%         .range, .resolution     - ADC specifications
%         .DCM_mounting           - Body to MAG mounting misalignment
%         .M_Deterministic        - Combined error matrix (Soft Iron + Nonorth + Mounting)
%         .biasHardIron           - Static hard iron bias [nT]
%         .sigmaWhite             - White noise RMS [nT]
%         .biasInstability        - Bias drift steady-state [nT]
%         .biasTau                - Bias correlation time [s]
%==========================================================================

    fprintf('\n=== Initializing MAG Model ===\n');
    
    % --- General Configuration ---
    MAG.rate = 10;           % Sampling frequency [Hz]
    MAG.dt   = 1 / MAG.rate; % Sampling time step [s]
    
    % Dynamic Range & Resolution       
    MAG.range      = 100000; % [nT] +/- full scale
    MAG.resolution = 0.5;    % [nT] ADC quantization step
    
    % Mounting misalignment errors (typical mechanical tolerances)
    rollErr  =  0.002; % Mounting error around X-axis [deg]
    pitchErr = -0.005; % Mounting error around Y-axis [deg]
    yawErr   =  0.003; % Mounting error around Z-axis [deg]
    
    % Construct misalignment DCM (small perturbation from ideal alignment)
    % This matrix is shared between gyro and accel (common mechanical mounting).
    MAG.DCM_mounting = angle2dcm(deg2rad(yawErr), ...
                                 deg2rad(pitchErr), ...
                                 deg2rad(rollErr), ...
                                 'ZYX');
    
    % --- Stochastic Noise Parameters ---
        MAG.noiseDensity = 0.1;               % White noise spectral density [nT/sqrt(Hz)]
        
        % Anti-alias / measurement bandwidth model for noise integration
        MAG.antiAlias.type = "onePole";       % "idealNyquist" | "onePole"
        MAG.antiAlias.fc   = 0.25 * MAG.rate; % [Hz]
        
        % Bias drift as bounded 1st-order Gauss-Markov process
        MAG.biasInstability = 5.0;            % Steady-state    [nT] (1σ, per axis)
        MAG.biasTau         = 3600;           % Correlation time [s]
    
    % --- Deterministic Errors ---
        MAG.hardIronLim      = 500.0;         % Static hard-iron magnitude limit [nT] (per axis, uniform)
        MAG.SF               = 1000;          % Diagonal scale errors [ppm] (1σ)
        MAG.softIronCoupling = 500;           % Off-diagonal symmetric coupling [ppm] (1σ)
        MAG.nonortho         = 0.1;           % Small non-orthogonality angles [deg] (1σ)
    
    % --- Derived Matrices & Pre-Allocation ---
    fprintf('\n=== Initializing MAG Model ===\n');
    
        % Soft iron / Scale factor (symmetric matrix near identity)
        SFErr  =               (MAG.SF * 1e-6) .* randn(3,1);
        SICErr = (MAG.softIronCoupling * 1e-6) .* randn(3,1); % xy, xz, yz couplings
        MAG.M_SoftIron = eye(3) + [ SFErr(1), SICErr(1), SICErr(2);
                                   SICErr(1),  SFErr(2), SICErr(3);
                                   SICErr(2), SICErr(3),  SFErr(3)];
         
        fprintf(' Soft-Iron / Scale-Factor Matrix: [%+.6f  %+.6f  %+.6f]\n', MAG.M_SoftIron(1,1), MAG.M_SoftIron(1,2), MAG.M_SoftIron(1,3));
        fprintf('                                  [%+.6f  %+.6f  %+.6f]\n', MAG.M_SoftIron(2,1), MAG.M_SoftIron(2,2), MAG.M_SoftIron(2,3));
        fprintf('                                  [%+.6f  %+.6f  %+.6f]\n', MAG.M_SoftIron(3,1), MAG.M_SoftIron(3,2), MAG.M_SoftIron(3,3));
        
        % Non-orthogonality (small-angle upper triangular correction)
        nonorthErr = deg2rad(MAG.nonortho) .* randn(3,1);
        MAG.M_Nonorth = [             1, -nonorthErr(3),  nonorthErr(2);
                         -nonorthErr(3),              1, -nonorthErr(1);
                          nonorthErr(2), -nonorthErr(1),              1];
        
        % Deterministic MAG matrix acting in MAG frame (after mounting)
        MAG.M_Deterministic = MAG.M_SoftIron * MAG.M_Nonorth * MAG.DCM_mounting;
    
    % Static hard-iron bias (MAG axes)
    MAG.biasHardIron = MAG.hardIronLim * (2*rand(3,1) - 1);
    
    fprintf(' Hard Iron Bias (MAG):            [%.2f, %.2f, %.2f] nT\n', MAG.biasHardIron(1), MAG.biasHardIron(2), MAG.biasHardIron(3));
    
    % White noise RMS per sample from density and assumed bandwidth
    if string(MAG.antiAlias.type) == "idealNyquist"
        BW = MAG.antiAlias.fc/2;
    elseif string(MAG.antiAlias.type) == "onePole"
        BW = 1.57 * MAG.antiAlias.fc; % ENBW for 1-pole low-pass
    else
        BW = MAG.antiAlias.fc/2;
    end
    MAG.sigmaWhite = MAG.noiseDensity * sqrt(BW);
    
    fprintf(' White Noise RMS:                 %.4f nT (per axis)\n', MAG.sigmaWhite);
    
    fprintf('=== MAG Initialization Complete ===\n');
end
