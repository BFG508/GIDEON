function STR = initializeSTR(nSTR)
%==========================================================================
% initializeSTR - Initialize Star Tracker hardware parameters
%
% Inputs:
%   nSTR - Number of star trackers (1 or 2)
%
% Outputs:
%   STR - Array of structures (1xnSTR) with fields:
%         .FOV                        - Field of view [deg]
%         .resolution                 - Sensor resolution [pixels]
%         .pixelSize                  - Pixel pitch [m]
%         .focalLength                - Optical focal length [m]
%         .centroidAccuracy           - Centroiding precision [pixels, 1σ]
%         .magnitudeLimit             - Visual magnitude threshold [-]
%         .rate                       - Measurement frequency [Hz]
%         .attitudeAccuracyCrossbore  - Cross-boresight accuracy [arcsec, 1σ]
%         .attitudeAccuracyRoll       - Roll-axis accuracy [arcsec, 1σ]
%         .maxAngularRate             - Maximum trackable rate [deg/s]
%         .DCM_B2S                    - Body to Sensor mounting DCM
%==========================================================================

    if nargin < 1
        nSTR = 2;  % Default: 2 star trackers
    end
    
    fprintf('\n=== Initializing Star Tracker Configuration (%d STRs) ===\n', nSTR);

    % Mounting misalignment errors (typical mechanical tolerances)
    rollErr  =  0.002; % Mounting error around X-axis [deg]
    pitchErr = -0.005; % Mounting error around Y-axis [deg]
    yawErr   =  0.003; % Mounting error around Z-axis [deg]
    
    % Construct misalignment DCM (small perturbation from ideal alignment)
    DCM_misalignment = angle2dcm(deg2rad(yawErr), ...
                                 deg2rad(pitchErr), ...
                                 deg2rad(rollErr), ...
                                 'ZYX');
    
    % Base STR specifications
    STR.FOV                       = 20;           % Field of view              [deg]
    STR.resolution                = [2048, 2048]; % Sensor resolution          [pixels]
    STR.pixelSize                 = 5.5e-6;       % Pixel pitch                [m]
    STR.focalLength               = 0.065;        % Optical focal length       [m]
    STR.centroidAccuracy          = 0.15;         % Centroiding precision      [pixels, 1σ]
    STR.magnitudeLimit            = 6.5;          % Absuolute magnitude threshold [-]
    STR.rate                      = 4;            % Measurement frequency      [Hz]
    STR.attitudeAccuracyCrossbore = 5;            % Cross-boresight accuracy   [arcsec, 1σ]
    STR.attitudeAccuracyRoll      = 20;           % Roll-axis accuracy         [arcsec, 1σ]
    STR.maxAngularRate            = 1.0;          % Maximum trackable rate     [deg/s]
    
    % STR1: Nominal boresight aligned with +Z body axis
    STR(1).DCM_B2S                   = DCM_misalignment;
    STR(1).FOV                       = STR.FOV;
    STR(1).resolution                = STR.resolution;
    STR(1).pixelSize                 = STR.pixelSize;
    STR(1).focalLength               = STR.focalLength;
    STR(1).centroidAccuracy          = STR.centroidAccuracy;
    STR(1).magnitudeLimit            = STR.magnitudeLimit;
    STR(1).rate                      = STR.rate;
    STR(1).attitudeAccuracyCrossbore = STR.attitudeAccuracyCrossbore;
    STR(1).attitudeAccuracyRoll      = STR.attitudeAccuracyRoll;
    
    fprintf(' STR1: Boresight +Z body (with %.3f° mounting error)\n', ...
                norm([rollErr, pitchErr, yawErr]));

    % STR2: 90° rotation around Y-axis for lateral coverage (if nSTR = 2)
    if nSTR == 2
        % Compose misalignment with nominal 90° Y-rotation
        DCM_STR2       = angle2dcm(-pi/2, 0, 0, 'YZX');
        STR(2)         = STR(1);
        STR(2).DCM_B2S = DCM_misalignment * DCM_STR2;

        fprintf(' STR2: Boresight +X body (90° Y-rotation + mounting error)\n');
    end

    fprintf('=== Star Tracker Initialization Complete ===\n');
end