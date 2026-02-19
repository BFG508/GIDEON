function fig = plotEKF_NavResults(t, truth, rEst, vEst, biasEst, PHist, imuMeas, saveFlag)
%==========================================================================
% plotEKF_NavResults: Comprehensive PV-EKF performance analysis and 
%                     visualization with translational errors, bias 
%                     estimation, uncertainty bounds, and consistency checks.
%
% INPUTS:
%   t        - Time vector                                           [s], 1xN
%   truth    - Ground truth structure (from generateGroundTruth) with fields:
%              .rECI       - True position in ECI                    [m], 3xN
%              .vECI       - True velocity in ECI                  [m/s], 3xN
%   rEst     - Estimated position history                            [m], 3xN
%   vEst     - Estimated velocity history                          [m/s], 3xN
%   biasEst  - Estimated accelerometer bias history                [m/s²], 3xN
%   PHist    - Covariance matrix history                                , 9x9xN
%   imuMeas  - IMU measurements structure (contains .accel.biasDyn)
%   saveFlag - (Optional) Boolean. If true, saves figures to 'Figures/EKF'.
%              Default: false
%
% OUTPUTS:
%   fig      - Array of figure handles for the generated plots (7x1)
%
% PLOTS GENERATED:
%   1) Position Error (x, y, z, Norm) with ±3σ bounds
%   2) Velocity Error (x, y, z, Norm) with ±3σ bounds
%   3) Accelerometer Bias estimation error with ±3σ bounds
%   4) Uncertainty evolution (3σ) for Pos, Vel, and Bias
%   5) NEES consistency check (chi-squared bounds)
%   6) 3D Orbital Trajectory (True vs Estimated)
%   7) Error statistics summary (bar charts)
%==========================================================================

    % Handle optional saveFlag
    if nargin < 7
        saveFlag = false;
    end
    
    nFig = 0;
    fig  = gobjects(7, 1);
    N    = length(t);
    tMin = t / 60;
    
    % Prepare output directory if saving
    if saveFlag
        saveDir = fullfile(pwd, 'Figures', 'EKF', 'Navigation');
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
            fprintf('\n✓ Created directory: %s\n', saveDir);
        end
        fprintf('\n--- Saving PV-EKF figures ---\n');
    end
    
    % Compute Errors
    rErr = rEst - truth.rECI;
    vErr = vEst - truth.vECI;
    
    % Bias error (Assuming baseline true dynamic bias is 0 for reference)
    biasErr  = (biasEst - imuMeas.accel.biasDyn) * 1000;
    
    % Extract uncertainty (1σ) from covariance
    sigmaPos  = zeros(3, N);
    sigmaVel  = zeros(3, N);
    sigmaBias = zeros(3, N);
    
    for k = 1:N
        sigmaPos(:,k)  = sqrt(diag(PHist(1:3, 1:3, k)));
        sigmaVel(:,k)  = sqrt(diag(PHist(4:6, 4:6, k)));
        sigmaBias(:,k) = sqrt(diag(PHist(7:9, 7:9, k))) * 1000;
    end
    
    % Steady state index for dynamic Y-axis scaling (ignore first 2%)
    idxConverged = max(1, round(0.02 * N)) : N;

    %% ------------------------------------------------------------------------
    % 1. POSITION ERROR
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'PV-EKF - Position Error', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [100, 100, 1000, 800]);

    tLayout = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    axisLabels = {'\Deltar_x [m]', '\Deltar_y [m]', '\Deltar_z [m]'};
    
    for i = 1:3
        ax = nexttile;
        hold(ax, 'on');
        
        plot(ax, tMin,        rErr(i,:),  'b-', 'LineWidth', 1.0, 'DisplayName', 'Error');
        plot(ax, tMin,  3*sigmaPos(i,:), 'r--', 'LineWidth', 1.2, 'DisplayName', '\pm3σ');
        plot(ax, tMin, -3*sigmaPos(i,:), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        
        ylabel(ax, axisLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
        set(ax, 'FontSize', 9);
        grid(ax, 'on');
        xlim(ax, [tMin(1), tMin(end)]);
        
        yMax = max(3 * sigmaPos(i, idxConverged));
        yMax = max(yMax * 1.5, 1); 
        ylim(ax, [-yMax, yMax]);
    end
    
    % Norm plot
    axNorm = nexttile;
    hold(axNorm, 'on');
    plot(axNorm, tMin, vecnorm(rErr, 2, 1), 'k-', 'LineWidth', 1.2, 'DisplayName', '||\Deltar||');

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
        saveFigure(fig(nFig), saveDir, 'EKF_Nav_fig1_position_error');
    end
    
    %% ------------------------------------------------------------------------
    % 2. VELOCITY ERROR
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'PV-EKF - Velocity Error', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [150, 100, 1000, 800]);

    tLayout = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    axisLabels = {'\Deltav_x [m/s]', '\Deltav_y [m/s]', '\Deltav_z [m/s]'};
    
    for i = 1:3
        ax = nexttile;
        hold(ax, 'on');
        
        plot(ax, tMin,        vErr(i,:),  'b-', 'LineWidth', 1.0, 'DisplayName', 'Error');
        plot(ax, tMin,  3*sigmaVel(i,:), 'r--', 'LineWidth', 1.2, 'DisplayName', '\pm3σ');
        plot(ax, tMin, -3*sigmaVel(i,:), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        
        ylabel(ax, axisLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
        set(ax, 'FontSize', 9);
        grid(ax, 'on');
        xlim(ax, [tMin(1), tMin(end)]);
        
        yMax = max(3 * sigmaVel(i, idxConverged));
        yMax = max(yMax * 1.5, 0.05); 
        ylim(ax, [-yMax, yMax]);
    end
    
    % Norm plot
    axNorm = nexttile;
    hold(axNorm, 'on');
    plot(axNorm, tMin, vecnorm(vErr, 2, 1), 'k-', 'LineWidth', 1.2, 'DisplayName', '||\Deltav||');

    ylabel(axNorm, '||\Deltav|| [m/s]', 'FontSize', 10, 'FontWeight', 'bold');
    grid(axNorm, 'on');
    set(axNorm, 'FontSize', 9);
    xlim(axNorm, [tMin(1), tMin(end)]);
    
    title(tLayout, 'Velocity Error', ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    lgd = legend(nexttile(1), 'Orientation', 'horizontal', 'FontSize', 11);
    lgd.Layout.Tile = 'north';
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'EKF_Nav_fig2_velocity_error');
    end
    
    %% ------------------------------------------------------------------------
    % 3. ACCELEROMETER BIAS ESTIMATION ERROR
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'PV-EKF - Accel Bias Error', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [200, 150, 900, 600]);

    tLayout = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    axisLabels = {'x [mm/s^2]', 'y [mm/s^2]', 'z [mm/s^2]'};
    
    for i = 1:3
        ax = nexttile;
        hold(ax, 'on');
        
        plot(ax, tMin,      biasErr(i,:),  'b-', 'LineWidth', 1.5, 'DisplayName', 'Error');
        plot(ax, tMin,  3*sigmaBias(i,:), 'r--', 'LineWidth', 1.2, 'DisplayName', '\pm3σ');
        plot(ax, tMin, -3*sigmaBias(i,:), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        
        ylabel(ax, axisLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
        set(ax, 'FontSize', 9);
        grid(ax, 'on');
        xlim(ax, [tMin(1), tMin(end)]);
        
        yMax = max(3 * sigmaBias(i, idxConverged));
        yMax = max(yMax * 1.5, 0.1); 
        ylim(ax, [-yMax, yMax]);
    end
    
    title(tLayout, 'Accelerometer Bias Error', ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    lgd = legend(nexttile(1), 'Orientation', 'horizontal', 'FontSize', 11);
    lgd.Layout.Tile = 'north';
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'EKF_Nav_fig3_accel_bias_error');
    end
    
    %% ------------------------------------------------------------------------
    % 4. POSITION & VELOCITY UNCERTAINTY EVOLUTION (3σ)
    % -------------------------------------------------------------------------
    
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'PV-EKF - Uncertainty Evolution', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [250, 200, 900, 800]);
                   
    tLayout = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % --- (1) Top: Position Uncertainty ---
    ax1 = nexttile;
    hold(ax1, 'on');
    plot(ax1, tMin, 3*sigmaPos(1,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'X');
    plot(ax1, tMin, 3*sigmaPos(2,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Y');
    plot(ax1, tMin, 3*sigmaPos(3,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Z');
    ylabel(ax1, '3\sigma Uncert. [m]', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax1, 'Position Uncertainty Evolution', 'FontSize', 12, 'FontWeight', 'bold');
    set(ax1, 'FontSize', 10); grid(ax1, 'on'); xlim(ax1, [tMin(1), tMin(end)]);
    ylim(ax1, [0, max(max(3*sigmaPos(:, idxConverged))) * 1.5]);
    legend(ax1, 'Location', 'northeast', 'FontSize', 10);
    
    % --- (2) Middle: Velocity Uncertainty ---
    ax2 = nexttile;
    hold(ax2, 'on');
    plot(ax2, tMin, 3*sigmaVel(1,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'X');
    plot(ax2, tMin, 3*sigmaVel(2,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Y');
    plot(ax2, tMin, 3*sigmaVel(3,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Z');
    ylabel(ax2, '3\sigma Uncert. [m/s]', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax2, 'Velocity Uncertainty Evolution', 'FontSize', 12, 'FontWeight', 'bold');
    set(ax2, 'FontSize', 10); grid(ax2, 'on'); xlim(ax2, [tMin(1), tMin(end)]);
    ylim(ax2, [0, max(max(3*sigmaVel(:, idxConverged))) * 1.5]);
    legend(ax2, 'Location', 'northeast', 'FontSize', 10);
    
    % --- (3) Bottom: Bias Uncertainty ---
    ax3 = nexttile;
    hold(ax3, 'on');
    plot(ax3, tMin, 3*sigmaBias(1,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'X');
    plot(ax3, tMin, 3*sigmaBias(2,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Y');
    plot(ax3, tMin, 3*sigmaBias(3,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Z');
    ylabel(ax3, '3\sigma Uncert. [mm/s^2]', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax3, 'Accel Bias Uncertainty Evolution', 'FontSize', 12, 'FontWeight', 'bold');
    set(ax3, 'FontSize', 10); grid(ax3, 'on'); xlim(ax3, [tMin(1), tMin(end)]);
    ylim(ax3, [0, max(max(3*sigmaBias(:, idxConverged))) * 1.5]);
    legend(ax3, 'Location', 'northeast', 'FontSize', 10);
    
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'EKF_Nav_fig4_uncertainty_evolution');
    end
    
    %% ------------------------------------------------------------------------
    % 5. NEES CONSISTENCY CHECK (6 DOF)
    % -------------------------------------------------------------------------
    
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'PV-EKF - NEES Consistency', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [500, 200, 1000, 700]);
    
    NEESnav = zeros(1, N);
    for k = 1:N
        % EKF translational error state (Position + Velocity)
        deltaX = [rErr(:,k); vErr(:,k)];
        
        % Covariance submatrix (Position + Velocity)
        Preg = PHist(1:6, 1:6, k) + eye(6)*1e-12;
        
        % NEES calculation
        NEESnav(k) = deltaX' * (Preg \ deltaX);
    end
    
    % Chi-squared bounds (95% confidence interval for 6 DOF)
    r1 = chi2inv(0.025, 6);
    r2 = chi2inv(0.975, 6);
    expectedNEES = 6;

    % Percentage inside bounds
    inside = sum(NEESnav > r1 & NEESnav < r2) / N * 100;

    tLayout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % --- (1) Top: Instantaneous NEES ---
    ax1 = nexttile;
    hold(ax1, 'on');
    plot(ax1, tMin, NEESnav, 'b-', 'LineWidth', 1.5, 'DisplayName', 'NEES');
    plot(ax1, tMin, r2*ones(size(t)), 'r--', 'LineWidth', 1.2, 'DisplayName', '95% CI');
    plot(ax1, tMin, r1*ones(size(t)), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    ylabel(ax1, 'NEES', 'FontSize', 10, 'FontWeight', 'bold');
    set(ax1, 'FontSize', 9); grid(ax1, 'on'); xlim(ax1, [tMin(1), tMin(end)]);
    
    % --- (2) Bottom: Cumulative Average NEES ---
    ax2 = nexttile;
    hold(ax2, 'on');
    avgNEES = cumsum(NEESnav) ./ (1:N);
    plot(ax2, tMin, avgNEES, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Average NEES');
    plot(ax2, tMin, expectedNEES*ones(size(t)), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Expected (6 DOF)');
    ylabel(ax2, 'Average NEES', 'FontSize', 10, 'FontWeight', 'bold');
    set(ax2, 'FontSize', 9); grid(ax2, 'on'); xlim(ax2, [tMin(1), tMin(end)]);
    
    title(tLayout, sprintf('Navigation Normalized Estimation Error Squared | %.1f%% inside 95%% CI', inside), ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'EKF_Nav_fig5_NEES_consistency');
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
    rNormMeas = vecnorm(rEst, 2, 1);
    latMeas   = asind(rEst(3,:) ./ rNormMeas);
    lonMeas   = atan2d(rEst(2,:), rEst(1,:));
    
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
    plot3(ax2,  rEst(1,:),  rEst(2,:),  rEst(3,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Measured');
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
        saveFigure(fig(nFig), saveDir, 'EKF_Nav_fig6_orbit');
    end

    %% ------------------------------------------------------------------------
    % 7. ERROR STATISTICS SUMMARY (BAR CHARTS)
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'PV-EKF - Error Statistics', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [400, 350, 1200, 600]);
    
    % Compute statistics
    stats.rMean = mean(abs(rErr), 2);
    stats.rStd  = std(rErr, 0, 2);
    stats.rMax  = max(abs(rErr), [], 2);
    stats.rRms  = sqrt(mean(rErr.^2, 2));
    
    stats.vMean = mean(abs(vErr), 2);
    stats.vStd  = std(vErr, 0, 2);
    stats.vMax  = max(abs(vErr), [], 2);
    stats.vRms  = sqrt(mean(vErr.^2, 2));
    
    stats.bMean = mean(abs(biasErr), 2);
    stats.bStd  = std(biasErr, 0, 2);
    stats.bMax  = max(abs(biasErr), [], 2);
    stats.bRms  = sqrt(mean(biasErr.^2, 2));
    
    x = 1:3; width = 0.2;
    
    % Position
    subplot(1,3,1); hold on; grid on;
    h(1) = bar(x - 1.5*width, stats.rMean, width, 'DisplayName', 'Mean');
    h(2) = bar(x - 0.5*width, stats.rStd, width, 'DisplayName', 'STD');
    h(3) = bar(x + 0.5*width, stats.rRms, width, 'DisplayName', 'RMS');
    h(4) = bar(x + 1.5*width, stats.rMax, width, 'DisplayName', 'Max');
    ylabel('Error [m]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Position Error', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'XTick', 1:3, 'XTickLabel', {'X', 'Y', 'Z'}, 'FontSize', 10);
    
    % Velocity
    subplot(1,3,2); hold on; grid on;
    bar(x - 1.5*width, stats.vMean, width);
    bar(x - 0.5*width, stats.vStd, width);
    bar(x + 0.5*width, stats.vRms, width);
    bar(x + 1.5*width, stats.vMax, width);
    ylabel('Error [m/s]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Velocity Error', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'XTick', 1:3, 'XTickLabel', {'X', 'Y', 'Z'}, 'FontSize', 10);
    
    % Accel Bias
    subplot(1,3,3); hold on; grid on;
    bar(x - 1.5*width, stats.bMean, width);
    bar(x - 0.5*width, stats.bStd, width);
    bar(x + 0.5*width, stats.bRms, width);
    bar(x + 1.5*width, stats.bMax, width);
    ylabel('Error [mm/s^2]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Accel Bias Error', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'XTick', 1:3, 'XTickLabel', {'X', 'Y', 'Z'}, 'FontSize', 10);
    
    % Shared Legend
    L = legend(h, 'Orientation', 'horizontal', 'FontSize', 11);
    L.Units = 'normalized';
    L.Position(1) = 0.5 - L.Position(3)/2;
    L.Position(2) = 0.95;
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'EKF_Nav_fig7_error_statistics');
    end
   
end