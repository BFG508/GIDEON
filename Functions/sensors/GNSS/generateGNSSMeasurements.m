function meas = generateGNSSMeasurements(t, rECI_true, vECI_true, qTrue, omegaTrue, GNSS)
%==========================================================================
% generateGNSSMeasurements: Generate synthetic GNSS PVT measurements
%                           (Position and Velocity) from truth data using 
%                           a space-grade receiver error model.
%
% DESCRIPTION:
%   Simulates a GNSS receiver output in the ECI frame. It models the
%   physical offset of the antenna (Lever Arm effect) coupled with attitude
%   dynamics, high-frequency receiver noise, and slowly varying correlated
%   errors (ephemeris, ionosphere, clock wander) using Gauss-Markov processes.
%
% INPUTS:
%   t               - Time vector                                       [s], 1xN
%   rECI_true       - True spacecraft CG position in ECI frame          [m], 3xN
%   vECI_true       - True spacecraft CG velocity in ECI frame        [m/s], 3xN
%   qTrue           - True attitude quaternion (ECI to Body)               , 4xN
%   omegaTrue       - True angular velocity in body frame           [rad/s], 3xN
%   GNSS            - GNSS parameter structure with required fields:
%                     .dt             - Sampling time                   [s]
%                     .leverArm       - Antenna offset from CG (Body)   [m], 3x1
%                     .sigmaPosWhite  - Position white noise 1σ         [m]
%                     .sigmaVelWhite  - Velocity white noise 1σ       [m/s]
%                     .tauPos         - Pos error correlation time      [s]
%                     .sigmaPosGM     - Pos Gauss-Markov steady-state   [m]
%                     .tauVel         - Vel error correlation time      [s]
%                     .sigmaVelGM     - Vel Gauss-Markov steady-state [m/s]
%
% OUTPUTS (struct meas):
%   meas.rECI       - Measured antenna position (Noisy)                 [m], 3xN
%   meas.vECI       - Measured antenna velocity (Noisy)               [m/s], 3xN
%   meas.posBiasDyn - Correlated position error (Ephemeris/Iono)        [m], 3xN
%   meas.velBiasDyn - Correlated velocity error (Clock drift)         [m/s], 3xN
%   meas.rClean     - True antenna position (CG + Lever Arm)            [m], 3xN
%   meas.vClean     - True antenna velocity (CG + Lever Arm)          [m/s], 3xN
%   meas.nSats      - Number of tracked satellites                      [-], 1xN
%
% ALGORITHM:
%   1. Lever Arm Kinematics:
%      rAnt_ECI = rCG_ECI + DCM_Body2ECI * LArm_body
%      vAnt_ECI = vCG_ECI + DCM_Body2ECI * (omegaBody x LArm_body)
%
%   2. Correlated Errors (1st-Order Gauss-Markov):
%      b_k+1 = exp(-dt/tau) * b_k + noiseDriving
%
%   3. Measurement Generation:
%      meas = cleanAntenna + correlatedBias + whiteNoise
%==========================================================================

    fprintf('\n=== Simulating GNSS Measurements ===\n');

    % --- Input Validation ---
    if nargin < 6
        error('generateGNSSMeasurements: Not enough input arguments.');
    end

    N = numel(t);
    if size(rECI_true, 2) ~= N || size(vECI_true, 2) ~= N
        error('Input truth arrays must have length N = length(t).');
    end

    % Pre-allocate output structure arrays
    meas.rECI       = zeros(3, N);
    meas.vECI       = zeros(3, N);
    meas.posBiasDyn = zeros(3, N);
    meas.velBiasDyn = zeros(3, N);
    meas.rClean     = zeros(3, N);
    meas.vClean     = zeros(3, N);
    meas.nSats      = zeros(1, N);

    % --- 1. Noise parameters & Gauss-Markov setup ---
    dt = GNSS.dt;

    % Position Gauss-Markov parameters
    if GNSS.tauPos <= 0
        aPos = 0; 
    else
        aPos = exp(-dt / GNSS.tauPos);
    end
    qPos = sqrt(max(1 - aPos^2, 0)) * GNSS.sigmaPosGM;
    posBias = zeros(3,1);

    % Velocity Gauss-Markov parameters
    if GNSS.tauVel <= 0
        aVel = 0; 
    else
        aVel = exp(-dt / GNSS.tauVel); 
    end
    qVel = sqrt(max(1 - aVel^2, 0)) * GNSS.sigmaVelGM;
    velBias = zeros(3,1);

    % --- 2. Main simulation loop ---
    for k = 1:N
        % 2.1 Kinematics & Lever arm effect
            % Get ECI to Body DCM and transpose to get Body to ECI
            DCM_ECI2B = quat2dcm(qTrue(:,k));
            DCM_B2ECI = DCM_ECI2B'; 
            
            % Antenna offset in ECI frame
            rLever_ECI = DCM_B2ECI * GNSS.leverArm;
            
            % Antenna velocity offset in ECI frame: v = w x r
            vLever_body = cross(omegaTrue(:,k), GNSS.leverArm);
            vLever_ECI  = DCM_B2ECI * vLever_body;
            
            % True antenna state
            rClean = rECI_true(:,k) + rLever_ECI;
            vClean = vECI_true(:,k) + vLever_ECI;
            
            meas.rClean(:,k) = rClean;
            meas.vClean(:,k) = vClean;

        % 2.2 Update correlated errors (Gauss-Markov)
        posBias = aPos * posBias + qPos * randn(3,1);
        velBias = aVel * velBias + qVel * randn(3,1);
        
        meas.posBiasDyn(:,k) = posBias;
        meas.velBiasDyn(:,k) = velBias;

        % 2.3 White noise generation
        posNoise = GNSS.sigmaPosWhite * randn(3,1);
        velNoise = GNSS.sigmaVelWhite * randn(3,1);

        % 2.4 Final PVT measurement
        meas.rECI(:,k) = rClean + posBias + posNoise;
        meas.vECI(:,k) = vClean + velBias + velNoise;
        
        % 2.5 Simulate number of tracked satellites (Realistic LEO: 6 to 12)
        % Slowly varying random walk rounded to nearest integer
        if k == 1
            nSats = 9.0;
        else
            nSats = nSats + 0.1 * randn();
            nSats = max(min(nSats, 12), 6); % Clamp between 6 and 12
        end
        meas.nSats(k) = round(nSats);
    end
    
    fprintf('  GNSS Lever Arm offset:  [%.2f, %.2f, %.2f] m\n', GNSS.leverArm);
    fprintf('  Average Pos Error:      %.2f m (White) + %.2f m (Correlated)\n', ...
            GNSS.sigmaPosWhite, GNSS.sigmaPosGM);
    fprintf('  Average Vel Error:      %.2f m/s (White) + %.2f m/s (Correlated)\n', ...
            GNSS.sigmaVelWhite, GNSS.sigmaVelGM);
    fprintf('  GNSS Data generation complete.\n');
end