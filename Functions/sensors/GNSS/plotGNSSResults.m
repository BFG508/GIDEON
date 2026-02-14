function fig = plotGNSSResults(t, truth, meas, GNSS, saveFlag)
%==========================================================================
% plotGNSSResults: Generate comprehensive visualizations of the GNSS PVT 
%                  measurements, including correlated noise analysis and 
%                  lever-arm kinematic coupling.
%
% INPUTS:
%   t        - Time vector                                      [s], 1xN
%   truth    - Ground truth structure containing .rECI, .vECI
%   meas     - GNSS measurements structure (from generateGNSSMeasurements)
%   GNSS     - GNSS parameters structure (for noise references)
%   saveFlag - (Optional) Boolean. If true, saves figures to 'Figures/GNSS'.
%              Default: false
%
% OUTPUTS:
%   fig      - Array of figure handles (5 figures).
%
% PLOTS GENERATED:
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
    t_min = t / 60;
    
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
    
    % Compute theoretical 3-sigma bounds for plots
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
        hold(ax, 'on'); grid(ax, 'on');
        
        plot(ax, t_min, errPos(i,:), 'b-', 'LineWidth', 1.0, 'DisplayName', 'Error');
        plot(ax, t_min,  sigmaPos3*ones(size(t)), 'r--', 'LineWidth', 1.2, 'DisplayName', '\pm3\sigma Expected');
        plot(ax, t_min, -sigmaPos3*ones(size(t)), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        
        ylabel(ax, axisLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
        set(ax, 'FontSize', 9);
        xlim(ax, [t_min(1), t_min(end)]);
        ylim(ax, [-sigmaPos3*1.5, sigmaPos3*1.5]);
    end
    
    % Norm Plot
    axNorm = nexttile;
    hold(axNorm, 'on'); grid(axNorm, 'on');
    plot(axNorm, t_min, errPosNorm, 'k-', 'LineWidth', 1.2, 'DisplayName', '||\Deltar||');
    ylabel(axNorm, '||\Deltar|| [m]', 'FontSize', 10, 'FontWeight', 'bold');
    set(axNorm, 'FontSize', 9);
    xlim(axNorm, [t_min(1), t_min(end)]);
    
    title(tLayout, 'GNSS Position Error (Measured vs True Antenna)', ...
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
        hold(ax, 'on'); grid(ax, 'on');
        
        plot(ax, t_min, errVel(i,:), 'b-', 'LineWidth', 1.0, 'DisplayName', 'Error');
        plot(ax, t_min,  sigmaVel3*ones(size(t)), 'r--', 'LineWidth', 1.2, 'DisplayName', '\pm3\sigma Expected');
        plot(ax, t_min, -sigmaVel3*ones(size(t)), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        
        ylabel(ax, axisLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
        set(ax, 'FontSize', 9);
        xlim(ax, [t_min(1), t_min(end)]);
        ylim(ax, [-sigmaVel3*1.5, sigmaVel3*1.5]);
    end
    
    % Norm Plot
    axNorm = nexttile;
    hold(axNorm, 'on'); grid(axNorm, 'on');
    plot(axNorm, t_min, errVelNorm, 'k-', 'LineWidth', 1.2, 'DisplayName', '||\Deltav||');
    ylabel(axNorm, '||\Deltav|| [m/s]', 'FontSize', 10, 'FontWeight', 'bold');
    set(axNorm, 'FontSize', 9);
    xlim(axNorm, [t_min(1), t_min(end)]);
    
    title(tLayout, 'GNSS Velocity Error (Measured vs True Antenna)', ...
          'FontSize', 14, 'FontWeight', 'bold');
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
    hold(ax1, 'on'); grid(ax1, 'on');
    plot(ax1, t_min, meas.posBiasDyn(1,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'X');
    plot(ax1, t_min, meas.posBiasDyn(2,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Y');
    plot(ax1, t_min, meas.posBiasDyn(3,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Z');
    ylabel(ax1, 'Pos Bias [m]', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax1, 'Correlated Position Error (e.g. Ephemeris / Ionosphere)', 'FontSize', 12, 'FontWeight', 'bold');
    set(ax1, 'FontSize', 9);
    xlim(ax1, [t_min(1), t_min(end)]);
    legend(ax1, 'Location', 'northeast', 'FontSize', 10);
    
    % Velocity Bias (Clock Wander)
    ax2 = nexttile;
    hold(ax2, 'on'); grid(ax2, 'on');
    plot(ax2, t_min, meas.velBiasDyn(1,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'X');
    plot(ax2, t_min, meas.velBiasDyn(2,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Y');
    plot(ax2, t_min, meas.velBiasDyn(3,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Z');
    ylabel(ax2, 'Vel Bias [m/s]', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax2, 'Correlated Velocity Error (e.g. Clock Wander)', 'FontSize', 12, 'FontWeight', 'bold');
    set(ax2, 'FontSize', 9);
    xlim(ax2, [t_min(1), t_min(end)]);
    legend(ax2, 'Location', 'northeast', 'FontSize', 10);
    
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'GNSS_fig3_gauss_markov_bias');
    end

    %% ------------------------------------------------------------------------
    % 4. LEVER ARM KINEMATIC DISTURBANCE
    % -------------------------------------------------------------------------
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'GNSS - Lever Arm Effect', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [250, 200, 900, 600]);
                   
    % Difference between CG and Antenna Phase Center
    leverPos = meas.rClean - truth.rECI;
    leverVel = meas.vClean - truth.vECI;
    
    tLayout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    ax1 = nexttile;
    hold(ax1, 'on'); grid(ax1, 'on');
    plot(ax1, t_min, leverPos(1,:), 'r-', 'LineWidth', 1.2, 'DisplayName', 'X');
    plot(ax1, t_min, leverPos(2,:), 'g-', 'LineWidth', 1.2, 'DisplayName', 'Y');
    plot(ax1, t_min, leverPos(3,:), 'b-', 'LineWidth', 1.2, 'DisplayName', 'Z');
    ylabel(ax1, 'Offset [m]', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax1, sprintf('Antenna Position Offset from Center of Mass (Lever Arm = %.2fm)', norm(GNSS.leverArm)), 'FontSize', 12, 'FontWeight', 'bold');
    set(ax1, 'FontSize', 9);
    xlim(ax1, [t_min(1), t_min(end)]);
    legend(ax1, 'Location', 'northeast', 'FontSize', 10);
    
    ax2 = nexttile;
    hold(ax2, 'on'); grid(ax2, 'on');
    plot(ax2, t_min, leverVel(1,:), 'r-', 'LineWidth', 1.2, 'DisplayName', 'X');
    plot(ax2, t_min, leverVel(2,:), 'g-', 'LineWidth', 1.2, 'DisplayName', 'Y');
    plot(ax2, t_min, leverVel(3,:), 'b-', 'LineWidth', 1.2, 'DisplayName', 'Z');
    ylabel(ax2, 'Vel Disturbance [m/s]', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax2, 'Induced Antenna Velocity due to Spacecraft Attitude Rates (\omega \times r)', 'FontSize', 12, 'FontWeight', 'bold');
    set(ax2, 'FontSize', 9);
    xlim(ax2, [t_min(1), t_min(end)]);
    legend(ax2, 'Location', 'northeast', 'FontSize', 10);
    
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'GNSS_fig4_lever_arm');
    end

    %% ------------------------------------------------------------------------
    % 5. CONSTELLATION TRACKING & 3D ORBIT
    % -------------------------------------------------------------------------
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'GNSS - Tracking & Orbit', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [300, 200, 1100, 500]);
                   
    % Tracking Plot
    subplot(1, 2, 1);
    stairs(t_min, meas.nSats, 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [min]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Number of Satellites', 'FontSize', 11, 'FontWeight', 'bold');
    title('Tracked GNSS Satellites (Visbility)', 'FontSize', 13, 'FontWeight', 'bold');
    ylim([0, 14]);
    xlim([t_min(1), t_min(end)]);
    set(gca, 'FontSize', 10);
    
    % 3D Orbit Plot
    subplot(1, 2, 2);
    REarth = 6378137; % [m]
    [x, y, z] = sphere(40);
    surf(x*REarth, y*REarth, z*REarth, 'FaceColor', [0.3 0.6 0.9], ...
         'EdgeColor', 'none', 'FaceAlpha', 0.2);
    hold on;
    
    plot3(truth.rECI(1,:), truth.rECI(2,:), truth.rECI(3,:), 'k-', 'LineWidth', 2, 'DisplayName', 'True Orbit');
    plot3(meas.rECI(1,:), meas.rECI(2,:), meas.rECI(3,:), 'r.', 'MarkerSize', 1, 'DisplayName', 'GNSS Meas');
    hold off;
    
    axis equal;
    grid on;
    xlabel('X [m]', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Y [m]', 'FontSize', 10, 'FontWeight', 'bold');
    zlabel('Z [m]', 'FontSize', 10, 'FontWeight', 'bold');
    title('Spacecraft Trajectory (ECI Frame)', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 10);
    view(45, 20);
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'GNSS_fig5_orbit_and_tracking');
    end

end