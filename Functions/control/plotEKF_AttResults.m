function fig = plotEKF_AttResults(t, truth, qEst, biasEst, PHist, imuMeas, saveFlag)
%==========================================================================
% plotEKF_AttResults: Generate quicklook figures for EKF truth vs measurements.
%
% INPUTS:
%   t        - Time vector                                           [s], 1xN
%   truth    - Ground truth structure (from generateGroundTruth) with fields:
%              .qTrue      - True attitude quaternion (ECI to Body)     , 4xN
%              .omegaTrue  - True angular velocity (body frame)  [rad/s], 3xN
%              .rECI       - True position in ECI                    [m], 3xN
%              .vECI       - True velocity in ECI                  [m/s], 3xN
%              .B_ECI      - True magnetic field in ECI             [nT], 3xN
%   qEst     - Estimated attitude quaternion history                    , 4xN
%   biasEst  - Estimated gyro bias history                       [rad/s], 3xN
%   PHist    - Covariance matrix history                                , 6x6xN
%   imuMeas  - IMU measurements structure (contains .gyro.omegaBody)
%   saveFlag - (Optional) Boolean. If true, saves figures to 'Figures/EKF'.
%              Default: false
%
% OUTPUTS:
%   fig      - Array of figure handles for the generated plots (8x1)
%
% PLOTS GENERATED:
%   1) Attitude error (3-axis Euler angles) with ±3σ bounds
%   2) Attitude error norm magnitude
%   3) Gyro bias estimation error with ±3σ bounds
%   4) Attitude and bias uncertainty evolution (3σ)
%   5) NEES consistency check (chi-squared bounds)
%   6) Quaternion components: true vs estimated
%   7) Angular velocity: true vs measured
%   8) Error statistics summary (bar charts)
%==========================================================================

    % Handle optional saveFlag
    if nargin < 7
        saveFlag = false;
    end
    
    nFig = 0;
    fig  = gobjects(8, 1);
    N    = length(t);
    
    % Prepare output directory if saving
    if saveFlag
        saveDir = fullfile(pwd, 'Figures', 'EKF', 'Attitude');
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
            fprintf('\n✓ Created directory: %s\n', saveDir);
        end
        fprintf('\n--- Saving EKF figures ---\n');
    end
    
    % Compute errors
    qError      = zeros(4, N);
    attErrorDeg = zeros(3, N);
    
    for k = 1:N
        % Quaternion error: qErr = qTrue^(-1) ⊗ qEst
        qError(:,k) = quatmultiply(quatinv(truth.qTrue(:,k)), qEst(:,k));
        
        % Convert to Euler angles (ZYX convention)
        DCM_error        = quat2dcm(qError(:,k));
        eulerErrorRad    = dcm2euler_ZYX(DCM_error);
        attErrorDeg(:,k) = rad2deg(eulerErrorRad);
    end
    
    % Gyro bias error (assume true bias = 0)
    biasTrue = zeros(3, N);
    biasErr  = rad2deg(biasEst - biasTrue) * 3600;
    
    % Extract uncertainty (1σ) from covariance
    sigmaAtt  = zeros(3, N);
    sigmaBias = zeros(3, N);
    
    for k = 1:N
        sigmaAtt(:,k)  = rad2deg(sqrt(diag(PHist(1:3, 1:3, k))));
        sigmaBias(:,k) = rad2deg(sqrt(diag(PHist(4:6, 4:6, k)))) * 3600;
    end
    
    %% ------------------------------------------------------------------------
    % 1. ATTITUDE ERROR (3-AXIS) WITH UNCERTAINTY BOUNDS
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'EKF - Attitude Error', ...
                       'Color', 'w', 'NumberTitle', 'off');

    tLayout = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    axisLabels = {'Roll [deg]', 'Pitch [deg]', 'Yaw [deg]'};
    
    idxConverged = max(1, round(0.02 * length(t))) : length(t);
    for i = 1:3
        ax = nexttile;
        hold(ax, 'on');
        
        % Plot data
        plot(ax, t/60, attErrorDeg(i,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Error');
        plot(ax, t/60,  3*sigmaAtt(i,:), 'r--', 'LineWidth', 1.2, 'DisplayName', '\pm3σ');
        plot(ax, t/60, -3*sigmaAtt(i,:), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        
        ylabel(ax, axisLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
        set(ax, 'FontSize', 9);
        grid(ax, 'on');
        xlim(ax, [t(1), t(end)]/60);
        
        yMax = max(3 * sigmaAtt(i, idxConverged));
        yMax = max(yMax * 1.5, 0.1); 
        ylim(ax, [-yMax, yMax]);
    end
    
    title(tLayout, 'Attitude Error (Euler Angles)', ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    % One legend for the whole figure
    lgd = legend(nexttile(1), 'Orientation', 'horizontal', 'FontSize', 11);
    lgd.Layout.Tile = 'north';
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'EKF_fig1_attitude_error');
    end
    
    %% ------------------------------------------------------------------------
    % 2. ATTITUDE ERROR NORM
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'EKF - Attitude Error Norm', ...
                       'Color', 'w', 'NumberTitle', 'off');
    
    attErrNorm = zeros(1, N);
    for k = 1:N
        attErrNorm(k) = 2 * acos(min(abs(qError(1,k)), 1));
    end
    meanErr = mean(rad2deg(attErrNorm));
    masErr  = max(rad2deg(attErrNorm));
    
    plot(t/60, rad2deg(attErrNorm), 'b-', 'LineWidth', 1.5);
    ylabel('Attitude Error [deg]', 'FontSize', 10, 'FontWeight', 'bold');
    title(sprintf('Attitude Error Magnitude | Mean = %.3f° | Max = %.3f°', meanErr, masErr), ...
          'FontSize', 13, 'FontWeight', 'bold');
    
    set(gca, 'FontSize', 10);
    grid on;
    xlim([t(1), t(end)]/60);
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'EKF_fig2_attitude_norm');
    end
    
    %% ------------------------------------------------------------------------
    % 3. GYRO BIAS ESTIMATION ERROR
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'EKF - Gyro Bias Error', ...
                       'Color', 'w', 'NumberTitle', 'off');

    tLayout = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    axisLabels = {'X-axis [deg/h]', 'Y-axis [deg/h]', 'Z-axis [deg/h]'};
    
    idxConverged = max(1, round(0.02 * length(t))) : length(t);
    for i = 1:3
        ax = nexttile;
        hold(ax, 'on');
        
        % Plot data
        plot(ax, t/60, biasErr(i,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Error');
        plot(ax, t/60,  3*sigmaBias(i,:), 'r--', 'LineWidth', 1.2, 'DisplayName', '\pm3σ');
        plot(ax, t/60, -3*sigmaBias(i,:), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        
        ylabel(ax, axisLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
        set(ax, 'FontSize', 9);
        grid(ax, 'on');
        xlim(ax, [t(1), t(end)]/60);
        
        yMax = max(3 * sigmaBias(i, idxConverged));
        yMax = max(yMax * 1.5, 0.1); 
        ylim(ax, [-yMax, yMax]);
    end
    
    title(tLayout, 'Gyroscope Bias Error', ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    % One legend for the whole figure
    lgd = legend(nexttile(1), 'Orientation', 'horizontal', 'FontSize', 11);
    lgd.Layout.Tile = 'north';
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'EKF_fig3_gyro_bias_error');
    end
    
    %% ------------------------------------------------------------------------
    % 4. ATTITUDE AND BIAS UNCERTAINTY EVOLUTION (3σ)
    % -------------------------------------------------------------------------
    
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'MEKF - Uncertainty Evolution', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [450, 450, 900, 600]);
                   
    % --- CREATE TILED LAYOUT ---
    tLayout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Define the steady-state index (ignore the first 2% of data for scaling)
    idxConverged = max(1, round(0.02 * length(t))) : length(t);
    
    % --- (1) Top: Attitude Uncertainty ---
    ax1 = nexttile;
    hold(ax1, 'on');
    
    plot(ax1, t/60, 3*sigmaAtt(1,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Roll');
    plot(ax1, t/60, 3*sigmaAtt(2,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Pitch');
    plot(ax1, t/60, 3*sigmaAtt(3,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Yaw');
    
    ylabel(ax1, '3σ Uncertainty [deg]', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax1, 'Attitude Uncertainty Evolution', 'FontSize', 12, 'FontWeight', 'bold');
    
    set(ax1, 'FontSize', 10); 
    grid(ax1, 'on');
    xlim(ax1, [t(1), t(end)]/60);

    yMax = max(3 * sigmaAtt(i, idxConverged));
    yMax = max(yMax * 1.5, 0.1); 
    ylim(ax1, [0, yMax]);
    
    legend(ax1, 'Location', 'northeast', 'FontSize', 10);
    
    % --- (2) Bottom: Bias Uncertainty ---
    ax2 = nexttile;
    hold(ax2, 'on');
    
    plot(ax2, t/60, 3*sigmaBias(1,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'X-axis');
    plot(ax2, t/60, 3*sigmaBias(2,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Y-axis');
    plot(ax2, t/60, 3*sigmaBias(3,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Z-axis');
    
    ylabel(ax2, '3σ Uncertainty [deg/h]', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax2, 'Gyroscope Bias Uncertainty Evolution', 'FontSize', 12, 'FontWeight', 'bold');
    
    set(ax2, 'FontSize', 10); 
    grid(ax2, 'on');
    xlim(ax2, [t(1), t(end)]/60);
    
    yMax = max(3 * sigmaBias(i, idxConverged));
    yMax = max(yMax * 1.5, 0.1); 
    ylim(ax2, [0, yMax]);
    
    legend(ax2, 'Location', 'northeast', 'FontSize', 10);
    
    % Global X Label
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'MEKF_fig4_uncertainty_evolution');
    end
    
    %% ------------------------------------------------------------------------
    % 5. NEES CONSISTENCY CHECK
    % -------------------------------------------------------------------------
    
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'EKF - NEES Consistency', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [500, 200, 1000, 700]);
    
    NEESatt = zeros(1, N);
    for k = 1:N
        % EKF error state from quaternion error
        qw         = max(abs(qError(1,k)), 1e-10);
        deltaTheta = 2 * qError(2:4, k) / qw;
        
        % Add regularization to avoid numerical issues
        Preg  = PHist(1:3, 1:3, k) + eye(3)*1e-12;
        
        % NEES = δθ' * P^(-1) * δθ
        NEESatt(k) = deltaTheta' * (Preg \ deltaTheta);
    end
    
    % Chi-squared bounds (95% confidence interval for 3 DOF)
    r1 = chi2inv(0.025, 3);
    r2 = chi2inv(0.975, 3);

    % Percentage inside bounds
    inside = sum(NEESatt > r1 & NEESatt < r2) / N * 100;

    tLayout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % --- (1) Top: Instantaneous NEES ---
    ax1 = nexttile;
    hold(ax1, 'on');
    
    plot(ax1, t/60, NEESatt, 'b-', 'LineWidth', 1.5, 'DisplayName', 'NEES');
    
    plot(ax1, t/60, r2*ones(size(t)), 'r--', 'LineWidth', 1.2, 'DisplayName', '95% CI');
    plot(ax1, t/60, r1*ones(size(t)), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    
    ylabel(ax1, 'NEES', 'FontSize', 10, 'FontWeight', 'bold');
    set(ax1, 'FontSize', 9);
    grid(ax1, 'on');
    xlim(ax1, [t(1), t(end)]/60);
    
    % --- (2) Bottom: Cumulative Average NEES ---
    ax2 = nexttile;
    hold(ax2, 'on');
    
    avgNEES = cumsum(NEESatt) ./ (1:N);
    plot(ax2, t/60, avgNEES, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Average NEES');
    plot(ax2, t/60, 3*ones(size(t)), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Expected (3 DOF)');
    
    ylabel(ax2, 'Average NEES', 'FontSize', 10, 'FontWeight', 'bold');
    set(ax2, 'FontSize', 9);
    grid(ax2, 'on');
    xlim(ax2, [t(1), t(end)]/60);
    
    % One legend for the whole figure
    title(tLayout, sprintf('Attitude Normalized Estimation Error Squared | %.1f%% inside 95%% CI', inside), ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    % lgd = legend(ax1, 'Orientation', 'horizontal', 'FontSize', 11);
    % lgd.Layout.Tile = 'north';
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'EKF_fig5_NEES_consistency');
    end
    
    %% ------------------------------------------------------------------------
    % 6. QUATERNION COMPONENTS: TRUE VS ESTIMATED
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'EKF - Quaternion Components', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [450, 50, 950, 950]);
    
    tLayout = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    labels = {'q_w', 'q_x', 'q_y', 'q_z'};
    
    for i = 1:4
        ax = nexttile;
        hold(ax, 'on'); grid(ax, 'on');
        
        % Plot data natively
        plot(ax, t/60, truth.qTrue(i,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'True');
        plot(ax, t/60, qEst(i,:), 'b--', 'LineWidth', 1.2, 'DisplayName', 'Estimated');
        
        ylabel(ax, labels{i}, 'FontSize', 10, 'FontWeight', 'bold');
        set(ax, 'FontSize', 9);
        xlim(ax, [t(1), t(end)]/60);
    end
    
    title(tLayout, 'Quaternion Components', ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    % One legend for the whole figure
    lgd = legend(nexttile(1), 'Orientation', 'horizontal', 'FontSize', 11);
    lgd.Layout.Tile = 'north';
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'EKF_fig6_quaternion_components');
    end
    
    %% ------------------------------------------------------------------------
    % 7. ANGULAR VELOCITY: TRUE VS ESTIMATED
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'EKF - Angular Velocity', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [450, 50, 950, 950]);
    
    tLayout = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    axisLabels = {'\omega_x [deg/s]', '\omega_y [deg/s]', '\omega_z [deg/s]'};
    
    for i = 1:3
        ax = nexttile;
        hold(ax, 'on'); 
        
        plot(ax, t/60, rad2deg(truth.omegaTrue(i,:)), 'r-', 'LineWidth', 1.5, ...
             'DisplayName', 'True');
        plot(ax, t/60, rad2deg(imuMeas.gyro.omegaBody(i,:)), 'b--', 'LineWidth', 1, ...
             'DisplayName', 'Measured');
        
        ylabel(ax, axisLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
        set(ax, 'FontSize', 9);
        grid(ax, 'on');
        xlim(ax, [t(1), t(end)]/60);
    end
    
    title(tLayout, 'Angular Velocity', ...
          'FontSize', 13, 'FontWeight', 'bold');
    xlabel(tLayout, 'Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    
    % One legend for the whole figure
    lgd = legend(nexttile(1), 'Orientation', 'horizontal', 'FontSize', 11);
    lgd.Layout.Tile = 'north';
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'EKF_fig7_angular_velocity');
    end
    
    %% ------------------------------------------------------------------------
    % 8. ERROR STATISTICS SUMMARY (BAR CHARTS)
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'EKF - Error Statistics', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [500, 150, 1000, 800]);
    
    % Compute statistics
    stats.attMean  = mean(abs(attErrorDeg), 2);
    stats.attStd   = std(attErrorDeg, 0, 2);
    stats.attMax   = max(abs(attErrorDeg), [], 2);
    stats.attRms   = sqrt(mean(attErrorDeg.^2, 2));
    
    stats.biasMean = mean(abs(biasErr), 2);
    stats.biasStd  = std(biasErr, 0, 2);
    stats.biasMax  = max(abs(biasErr), [], 2);
    stats.biasRms  = sqrt(mean(biasErr.^2, 2));
    
    % Plot bar charts
    subplot(1,2,1);
    x     = 1:3;
    width = 0.2;
    
    h(1) = bar(x - 1.5*width, stats.attMean, width, 'DisplayName', 'Mean');
    hold on;
    h(2) = bar(x - 0.5*width, stats.attStd, width, 'DisplayName', 'STD');
    h(3) = bar(x + 0.5*width, stats.attRms, width, 'DisplayName', 'RMS');
    h(4) = bar(x + 1.5*width, stats.attMax, width, 'DisplayName', 'Max');
    
    ylabel('Error [deg]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Attitude Error Statistics', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    set(gca, 'XTick', 1:3, 'XTickLabel', {'Roll', 'Pitch', 'Yaw'}, 'FontSize', 10);
    
    subplot(1,2,2);
    bar(x - 1.5*width, stats.biasMean, width);
    hold on;
    bar(x - 0.5*width, stats.biasStd, width);
    bar(x + 0.5*width, stats.biasRms, width);
    bar(x + 1.5*width, stats.biasMax, width);
    
    ylabel('Error [deg/h]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Gyro Bias Error Statistics', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    set(gca, 'XTick', 1:3, 'XTickLabel', {'X', 'Y', 'Z'}, 'FontSize', 10);
    
    % One legend for the whole figure
    L = legend(h, 'Orientation', 'horizontal', 'FontSize', 11);
    L.Units = 'normalized';
    L.Position(1) = 0.5 - L.Position(3)/2;
    L.Position(2) = 0.955;
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'EKF_fig8_error_statistics');
    end
   
end
