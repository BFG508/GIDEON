function truth = generateGroundTruth(t, epoch, orbitalElements, scParams, attParams)
%==========================================================================
% generateTruthTrajectory - Generate realistic spacecraft trajectory by
%                           orchestrating orbital and attitude simulations.
%
% Inputs:
%   t               - Time vector                          [s], 1xN
%   epoch           - Simulation UTC start time            (datetime object)
%   orbitalElements - Struct with initial orbital parameters:
%       .SMA        - Semi-major axis                      [m]
%       .ECC        - Eccentricity [0, 1)                  [-]
%       .INC        - Inclination                        [deg]
%       .RAAN       - Right Ascension of the Asc. Node   [deg]
%       .AOP        - Argument of Perigee                [deg]
%       .TA         - True Anomaly at t = 0              [deg]
%       .epoch      - Datetime object (UTC)
%   scParams        - Struct with spacecraft physical parameters:
%       .mass       - Spacecraft total mass               [kg]
%       .areaDrag   - Cross-sectional area for drag       [m²]
%       .areaSRP    - Cross-sectional area for SRP        [m²]
%       .Cd         - Atmospheric drag coefficient         [-]
%       .Cr         - Solar radiation pressure coeff.      [-]
%   attParams       - Struct with attitude control parameters:
%       .I          - Inertia tensor (3x3)             [kg·m²]
%       .dipole     - Magnetic dipole moment (3x1)      [A·m²]
%       .Kp         - Proportional gain (PD control)     [N·m]
%       .Kd         - Derivative gain (PD control)     [N·m·s]
%       .mode       - Pointing mode (e.g., 'LVLH')
%
% Outputs:
%   truth - Structure with fields:
%       .t          - Time vector                            [s], 1xN
%       .rECI       - Position in ECI                        [m], 3xN
%       .vECI       - Velocity in ECI                      [m/s], 3xN
%       .qTrue      - Attitude quaternion (ECI to Body)      [-], 4xN
%       .omegaTrue  - Angular velocity (Body frame)      [rad/s], 3xN
%       .B_ECI      - Magnetic field in ECI                 [nT], 3xN
%       .torques    - Applied torques (Body frame)         [N·m], struct
%==========================================================================

    fprintf('\n=== Generating Truth Trajectory ===\n');
    t = t(:)'; % Ensure row vector
    N = numel(t);
    
    %% 1. Orbital propagation
    fprintf(' Propagating orbit ... ');
    [rECI, vECI] = simulateOrbit(t, orbitalElements, scParams);
    fprintf('Done.\n');
    
    %% 2. Geomagnetic Field (IGRF-13)
    fprintf(' Simulating geomagnetic field (IGRF-13) ... ');
    B_ECI = simulateGeomagneticField(t, rECI, epoch);
    fprintf('Done.\n');
    
    %% 3. Sun Ephemeris
    fprintf(' Computing Sun ephemeris ... ');
    sunPos = zeros(3,N);
    for k = 1:N
        sunPos(:,k) = sunPosition(epoch + seconds(t(k)));
    end
    fprintf('Done.\n');
    
    %% 4. Attitude Dynamics
    fprintf(' Integrating attitude dynamics ... ');
    [qTrue, omegaTrue, torques] = simulateAttitude(t, rECI, vECI, B_ECI, sunPos, attParams, scParams, epoch);
    fprintf('Done.\n');
    
    %% 5. Assemble Output
    truth.t         = t;
    truth.rECI      = rECI;
    truth.vECI      = vECI;
    truth.qTrue     = qTrue;
    truth.omegaTrue = omegaTrue;
    truth.B_ECI     = B_ECI;
    truth.torques   = torques;
    
    fprintf('=== Truth Trajectory Complete ===\n');
end