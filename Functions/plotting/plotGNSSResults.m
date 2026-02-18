function fig = plotGNSSResults(t, truth, meas, GNSS, saveFlag)
%==========================================================================
% plotGNSSResults: Generate comprehensive visualizations of the GNSS PVT 
%                  measurements, including correlated noise analysis and 
%                  lever-arm kinematic coupling.
%
% Inputs:
%   t        - Time vector                                                 [s], 1xN
%   truth    - Ground truth structure containing:
%              * .rECI          : True center of mass position (ECI)       [m], 3xN
%              * .vECI          : True center of mass velocity (ECI)     [m/s], 3xN
%   meas     - GNSS measurements structure containing:
%              * .rECI          : Measured antenna position (ECI)           [m], 3xN
%              * .rClean        : True antenna phase center pos (ECI)       [m], 3xN
%              * .vECI          : Measured antenna velocity (ECI)         [m/s], 3xN
%              * .vClean        : True antenna phase center vel (ECI)     [m/s], 3xN
%              * .posBiasDyn    : Gauss-Markov correlated pos bias          [m], 3xN
%              * .velBiasDyn    : Gauss-Markov correlated vel bias        [m/s], 3xN
%              * .nSats         : Number of visible/tracked satellites         , 1xN
%   GNSS     - GNSS configuration parameters structure containing:
%              * .sigmaPosWhite : White noise std deviation for position    [m]
%              * .sigmaPosGM    : Gauss-Markov std deviation for position   [m]
%              * .sigmaVelWhite : White noise std deviation for velocity  [m/s]
%              * .sigmaVelGM    : Gauss-Markov std deviation for velocity [m/s]
%              * .leverArm      : Antenna offset from Center of Mass        [m], 3x1
%   saveFlag - (Optional) Boolean. If true, saves figures to 'Figures/GNSS'.
%              Default: false
%
% Outputs:
%   fig      - Array of figure handles (5x1)
%
% Plots generated:
%   1) Position Error (X, Y, Z, Norm) vs Antenna Truth
%   2) Velocity Error (X, Y, Z, Norm) vs Antenna Truth
%   3) Gauss-Markov Correlated Bias Evolution
%   4) Lever Arm Kinematic Disturbance (Attitude-Translation coupling)
%   5) Visible Satellites and 3D Orbit trajectory
%==========================================================================

    % Handle optional saveFlag
    if nargin < 5
        saveFlag = false;
    end
    
    nFig = 0;
    fig  = gobjects(5, 1);
    
    % Time in minutes for better X-axis readability
    tMin = t / 60;
    
    % Prepare output directory if saving
    if saveFlag
        saveDir = fullfile(pwd, 'Figures', 'GNSS');
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
            fprintf('\n✓ Created directory: %s\n', saveDir);
        end
        fprintf('\n--- Saving GNSS figures ---\n');
    end
    
    % Compute pure sensor errors (Measured vs True Antenna Phase Center)
    errPos     = meas.rECI - meas.rClean;
    errPosNorm = vecnorm(errPos, 2, 1);
    
    errVel     = meas.vECI - meas.vClean;
    errVelNorm = vecnorm(errVel, 2, 1);
    
    % Compute theoretical 3σ bounds for plots
    sigmaPos3 = 3 * sqrt(GNSS.sigmaPosWhite^2 + GNSS.sigmaPosGM^2);
    sigmaVel3 = 3 * sqrt(GNSS.sigmaVelWhite^2 + GNSS.sigmaVelGM^2);

    %% ------------------------------------------------------------------------
    % 1. POSITION MEASUREMENT ERROR
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'GNSS - Position Error', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [100, 100, 1000, 800]);
                   
    tLayout = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    axisLabels = {'\Deltar_x [m]', '\Deltar_y [m]', '\Deltar_z [m]'};
    
    for i = 1:3
        ax = nexttile;
        hold(ax, 'on');
        
        plot(ax, tMin,              errPos(i,:), 'b-', 'LineWidth', 1.0, 'DisplayName', 'Error');
        plot(ax, tMin,  sigmaPos3*ones(size(t)), 'r--', 'LineWidth', 1.2, 'DisplayName', '\pm3σ');
        plot(ax, tMin, -sigmaPos3*ones(size(t)), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        
        ylabel(ax, axisLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
        set(ax, 'FontSize', 9);
         grid(ax, 'on');
        xlim(ax, [tMin(1), tMin(end)]);
        ylim(ax, [-sigmaPos3*1.5, sigmaPos3*1.5]);
    end
    
    % Norm plot
    axNorm = nexttile;
    hold(axNorm, 'on');
    plot(axNorm, tMin, errPosNorm, 'k-', 'LineWidth', 1.2, 'DisplayName', '||\Deltar||');
    ylabel(axNorm, '||\Deltar|| [m]', 'FontSize', 10, 'FontWeight', 'bold');

    set(axNorm, 'FontSize', 9);
    grid(axNorm, 'on');
    xlim(axNorm, [tMin(1), tMin(end)]);
    
    title(tLayout, 'Position Error', ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    lgd = legend(nexttile(1), 'Orientation', 'horizontal', 'FontSize', 11);
    lgd.Layout.Tile = 'north';
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'GNSS_fig1_position_error');
    end

    %% ------------------------------------------------------------------------
    % 2. VELOCITY MEASUREMENT ERROR
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'GNSS - Velocity Error', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [150, 100, 1000, 800]);
                   
    tLayout = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    axisLabels = {'\Deltav_x [m/s]', '\Deltav_y [m/s]', '\Deltav_z [m/s]'};
    
    for i = 1:3
        ax = nexttile;
        hold(ax, 'on');
        
        plot(ax, tMin,              errVel(i,:), 'b-', 'LineWidth', 1.0, 'DisplayName', 'Error');
        plot(ax, tMin,  sigmaVel3*ones(size(t)), 'r--', 'LineWidth', 1.2, 'DisplayName', '\pm3σ');
        plot(ax, tMin, -sigmaVel3*ones(size(t)), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        
        ylabel(ax, axisLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
        set(ax, 'FontSize', 9);
        grid(ax, 'on');
        xlim(ax, [tMin(1), tMin(end)]);
        ylim(ax, [-sigmaVel3*1.5, sigmaVel3*1.5]);
    end
    
    % Norm Plot
    axNorm = nexttile;
    hold(axNorm, 'on');
    plot(axNorm, tMin, errVelNorm, 'k-', 'LineWidth', 1.2, 'DisplayName', '||\Deltav||');
    ylabel(axNorm, '||\Deltav|| [m/s]', 'FontSize', 10, 'FontWeight', 'bold');

    set(axNorm, 'FontSize', 9); 
    grid(axNorm, 'on');
    xlim(axNorm, [tMin(1), tMin(end)]);
    
    title(tLayout, 'Velocity Error', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    lgd = legend(nexttile(1), 'Orientation', 'horizontal', 'FontSize', 11);
    lgd.Layout.Tile = 'north';
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'GNSS_fig2_velocity_error');
    end

    %% ------------------------------------------------------------------------
    % 3. GAUSS-MARKOV CORRELATED BIAS EVOLUTION
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'GNSS - Correlated Bias (Gauss-Markov)', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [200, 200, 900, 600]);
                   
    tLayout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Position Bias (Ephemeris, Ionosphere)
    ax1 = nexttile;
    hold(ax1, 'on');
    plot(ax1, tMin, meas.posBiasDyn(1,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'x');
    plot(ax1, tMin, meas.posBiasDyn(2,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'y');
    plot(ax1, tMin, meas.posBiasDyn(3,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'z');
    ylabel(ax1, 'Position Bias [m]', 'FontSize', 10, 'FontWeight', 'bold');
    
    grid(ax1, 'on');
    set(ax1, 'FontSize', 9);
    xlim(ax1, [tMin(1), tMin(end)]);
    
    % Velocity Bias (Clock Wander)
    ax2 = nexttile;
    hold(ax2, 'on');
    plot(ax2, tMin, meas.velBiasDyn(1,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'X');
    plot(ax2, tMin, meas.velBiasDyn(2,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Y');
    plot(ax2, tMin, meas.velBiasDyn(3,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Z');
    ylabel(ax2, 'Velocity Bias [m/s]', 'FontSize', 10, 'FontWeight', 'bold');
    set(ax2, 'FontSize', 9);
    grid(ax2, 'on');
    xlim(ax2, [tMin(1), tMin(end)]);
    
    title(tLayout, 'Correlated Error', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    lgd = legend(nexttile(1), 'Orientation', 'horizontal', 'FontSize', 11);
    lgd.Layout.Tile = 'north';
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'GNSS_fig3_gauss_markov_bias');
    end

    %% ------------------------------------------------------------------------
    % 4. LEVER ARM KINEMATIC DISTURBANCE
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'GNSS - Lever Arm Effect', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [200, 200, 900, 600]);
                   
    % Difference between CG and Antenna Phase Center
    leverPos = meas.rClean - truth.rECI;
    leverVel = meas.vClean - truth.vECI;
    
    tLayout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    ax1 = nexttile;
    hold(ax1, 'on');
    plot(ax1, tMin, leverPos(1,:), 'r-', 'LineWidth', 1.2, 'DisplayName', 'x');
    plot(ax1, tMin, leverPos(2,:), 'g-', 'LineWidth', 1.2, 'DisplayName', 'y');
    plot(ax1, tMin, leverPos(3,:), 'b-', 'LineWidth', 1.2, 'DisplayName', 'z');

    ylabel(ax1, 'Position Offset [m]', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax1, sprintf('Antenna Position Offset from Center of Mass (Lever Arm = %.2fm)', ...
          norm(GNSS.leverArm)), 'FontSize', 12, 'FontWeight', 'bold');
    set(ax1, 'FontSize', 9);
    grid(ax1, 'on');
    xlim(ax1, [tMin(1), tMin(end)]);
    
    ax2 = nexttile;
    hold(ax2, 'on');
    plot(ax2, tMin, leverVel(1,:), 'r-', 'LineWidth', 1.2, 'DisplayName', 'X');
    plot(ax2, tMin, leverVel(2,:), 'g-', 'LineWidth', 1.2, 'DisplayName', 'Y');
    plot(ax2, tMin, leverVel(3,:), 'b-', 'LineWidth', 1.2, 'DisplayName', 'Z');
    ylabel(ax2, 'Velocity Disturbance [m/s]', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax2, 'Induced Antenna Velocity due to Spacecraft Attitude Rates', ...
          'FontSize', 12, 'FontWeight', 'bold');
    set(ax2, 'FontSize', 9);
    grid(ax2, 'on');
    xlim(ax2, [tMin(1), tMin(end)]);
    
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    lgd = legend(ax1, 'Orientation', 'horizontal', 'FontSize', 11);
    lgd.Layout.Tile = 'north';
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'GNSS_fig4_lever_arm');
    end

    %% ------------------------------------------------------------------------
    % 5. CONSTELLATION TRACKING
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'GNSS - Constellation Tracking', ...
                       'Color', 'w', 'NumberTitle', 'off');
                   
    stairs(tMin, meas.nSats, 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [min]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Number of Satellites', 'FontSize', 11, 'FontWeight', 'bold');
    title('Tracked GNSS Satellites (Visibility)', 'FontSize', 13, 'FontWeight', 'bold');
    ylim([0, 14]);
    xlim([tMin(1), tMin(end)]);
    set(gca, 'FontSize', 10);
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'GNSS_fig5_tracking');
    end

    %% ------------------------------------------------------------------------
    % 6. 2D GROUND TRACK & 3D ORBIT
    % -------------------------------------------------------------------------
    
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'GNSS - 2D Projection & 3D Orbit', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [100, 200, 1400, 600]);
                   
    REarth = 6378137.0; % [m]
    
    % Try to load the image once to use in both plots
    hasImage = false;
    try
        earthImg = imread('Data\blueMarble.jpg');
        earthImg = imresize(earthImg, [1024, 2048]);
        hasImage = true;
    catch
        warning('Could not find "blueMarble.jpg".');
    end
    
    % ---------------------------------------------------------
    % Subplot 1: 2D Projection
    % ---------------------------------------------------------
    ax1 = nexttile;
    hold(ax1, 'on');
    
    % Draw 2D map if image is available
    if hasImage
        % Map image to Longitude [-180, 180] and Latitude [-90, 90]
        imagesc(ax1, [-180 180], [-90 90], flipud(earthImg));
        set(ax1, 'YDir', 'normal'); % Ensure North is at the top
    end
    
    % ECI to Spherical (RA / Dec) transformation for ground track
    % True Orbit
    rNormTrue = vecnorm(truth.rECI, 2, 1);
    latTrue   = asind(truth.rECI(3,:) ./ rNormTrue);
    lonTrue   = atan2d(truth.rECI(2,:), truth.rECI(1,:));
    
    % GNSS Measurements
    rNormMeas = vecnorm(meas.rECI, 2, 1);
    latMeas   = asind(meas.rECI(3,:) ./ rNormMeas);
    lonMeas   = atan2d(meas.rECI(2,:), meas.rECI(1,:));
    
    % Insert NaNs to prevent horizontal lines when crossing the +/- 180 deg boundary
    wrapIdxTrue = find(abs(diff(lonTrue)) > 180);
    lonTrue(wrapIdxTrue) = NaN;
    latTrue(wrapIdxTrue) = NaN;
    
    wrapIdxMeas = find(abs(diff(lonMeas)) > 180);
    lonMeas(wrapIdxMeas) = NaN;
    latMeas(wrapIdxMeas) = NaN;
    
    plot(ax1, lonTrue, latTrue, 'g-', 'LineWidth', 2, 'DisplayName', 'Truth');
    plot(ax1, lonMeas, latMeas, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Measured');
    
    axis(ax1, 'equal');
    axis(ax1, 'tight');
    xlim(ax1, [-180 180]);
    ylim(ax1, [-90 90]);
    grid(ax1, 'on');
    xlabel(ax1, 'Longitude [deg]', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel(ax1, 'Latitude [deg]', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax1, '2D Trajectory Projection', 'FontSize', 13, 'FontWeight', 'bold');
    
    % ---------------------------------------------------------
    % Subplot 2: 3D ORBIT
    % ---------------------------------------------------------
    ax2 = nexttile;
    [x, y, z] = sphere(60);
    hold(ax2, 'on');
    if hasImage
        % Flip the image vertically to align with the 3D Z-axis
        surf(ax2, x*REarth, y*REarth, z*REarth, ...
             'CData', flipud(earthImg), ...
             'FaceColor', 'texturemap', ...
             'EdgeColor', 'none', ...
             'FaceLighting', 'gouraud', ...
             'AmbientStrength', 0.2, ...   
             'DiffuseStrength', 0.9, ...
             'HandleVisibility', 'off');
    else
        surf(ax2, x*REarth, y*REarth, z*REarth, 'FaceColor', [0.3 0.6 0.9], ...
             'EdgeColor', 'none', 'FaceAlpha', 0.2, ...
             'HandleVisibility', 'off');
    end
    plot3(ax2, truth.rECI(1,:), truth.rECI(2,:), truth.rECI(3,:), 'g-', 'LineWidth', 2, 'DisplayName', 'Truth');
    plot3(ax2,  meas.rECI(1,:),  meas.rECI(2,:),  meas.rECI(3,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Measured');
    hold(ax2, 'off');
    
    axis(ax2, 'equal');
    xlabel(ax2, 'x [m]', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel(ax2, 'y [m]', 'FontSize', 10, 'FontWeight', 'bold');
    zlabel(ax2, 'z [m]', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax2, '3D Spacecraft Trajectory', 'FontSize', 13, 'FontWeight', 'bold');
    grid(ax2, 'on');
    view(ax2, 45, 20);
    
    lgd = legend(ax1, 'Orientation', 'horizontal', 'FontSize', 11);
    lgd.Layout.Tile = 'north';
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'GNSS_fig6_orbit');
    end

end