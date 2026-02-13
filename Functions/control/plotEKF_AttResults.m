function fig = plotEKF_AttResults(t, truth, qEst, biasEst, P_hist, imuMeas, saveFlag)
%==========================================================================
% plotEKF_AttResults: Comprehensive MEKF performance analysis and visualization
%                     with attitude errors, bias estimation, uncertainty bounds,
%                     and consistency checks.
%
% Inputs:
%   t        - Time vector                                        [s], 1xN
%   truth    - Truth structure with fields:
%              .qTrue      - True quaternion                          , 4xN
%              .omegaTrue  - True angular velocity            [rad/s], 3xN
%              .rECI, .vECI, .B_ECI, etc.
%   qEst     - Estimated quaternion                                   , 4xN
%   biasEst  - Estimated gyro bias                            [rad/s], 3xN
%   P_hist   - Covariance history                                 , 6x6xN
%   imuMeas  - IMU measurements structure
%   saveFlag - (Optional) Boolean. If true, saves figures to 'Figures/MEKF'.
%              Default: false
%
% Outputs:
%   fig      - Array of figure handles
%
% Plots generated:
%   1) Attitude error (3-axis Euler angles) with ±3σ bounds
%   2) Attitude error norm and quaternion normalization check
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
        saveDir = fullfile(pwd, 'Figures', 'MEKF');
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
            fprintf('\n✓ Created directory: %s\n', saveDir);
        end
        fprintf('\n--- Saving MEKF figures ---\n');
    end
    
    %% ------------------------------------------------------------------------
    % COMPUTE ERRORS
    % -------------------------------------------------------------------------
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
    biasTrue       = zeros(3, N);
    biasErrorDegph = rad2deg(biasEst - biasTrue) * 3600;
    
    % Extract uncertainty (1σ) from covariance
    sigmaAtt  = zeros(3, N);
    sigmaBias = zeros(3, N);
    
    for k = 1:N
        sigmaAtt(:,k)  = rad2deg(sqrt(diag(P_hist(1:3, 1:3, k))));
        sigmaBias(:,k) = rad2deg(sqrt(diag(P_hist(4:6, 4:6, k)))) * 3600;
    end
    
    %% ------------------------------------------------------------------------
    % 1. ATTITUDE ERROR (3-AXIS) WITH UNCERTAINTY BOUNDS
    % -------------------------------------------------------------------------
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'MEKF - Attitude Error', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [100, 100, 1370, 890]);
    
    axisLabels = {'Roll (X) [deg]', 'Pitch (Y) [deg]', 'Yaw (Z) [deg]'};
    
    for i = 1:3
        subplot(3,1,i);
        hold on; grid on;
        
        plot(t/60, attErrorDeg(i,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Error');
        plot(t/60,  3*sigmaAtt(i,:), 'r--', 'LineWidth', 1.2, 'DisplayName', '±3σ');
        plot(t/60, -3*sigmaAtt(i,:), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        
        ylabel(axisLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 9);
        
        if i == 1
            title('Attitude Estimation Error (Euler Angles)', ...
                  'FontSize', 13, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
        end
        
        set(gca, 'FontSize', 9);
        xlim([t(1)/60, t(end)/60]);
        
        rmsError = rms(attErrorDeg(i,:));
        text(0.02, 0.95, sprintf('RMS = %.3f°', rmsError), ...
             'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'white');
    end
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'MEKF_fig1_attitude_error');
    end
    
    %% ------------------------------------------------------------------------
    % 2. ATTITUDE ERROR NORM & QUATERNION NORMALIZATION CHECK
    % -------------------------------------------------------------------------
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'MEKF - Attitude Error Norm', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [150, 150, 1370, 600]);
    
    subplot(2,1,1);
    hold on; grid on;
    
    attErrorNorm = zeros(1, N);
    for k = 1:N
        attErrorNorm(k) = 2 * acos(min(abs(qError(1,k)), 1));
    end
    
    plot(t/60, rad2deg(attErrorNorm), 'b-', 'LineWidth', 1.5);
    ylabel('Attitude Error [deg]', 'FontSize', 10, 'FontWeight', 'bold');
    title('Total Attitude Error Magnitude', 'FontSize', 13, 'FontWeight', 'bold');
    
    set(gca, 'FontSize', 9);
    xlim([t(1)/60, t(end)/60]);
    
    meanErr = mean(rad2deg(attErrorNorm));
    maxErr  = max(rad2deg(attErrorNorm));
    text(0.02, 0.95, sprintf('Mean = %.3f°\nMax = %.3f°', meanErr, maxErr), ...
         'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'white', ...
         'VerticalAlignment', 'top');
    
    subplot(2,1,2);
    hold on; grid on;
    
    qNormEst  = vecnorm(qEst, 2, 1);
    qNormTrue = vecnorm(truth.qTrue, 2, 1);
    
    plot(t/60, qNormEst - 1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Estimated');
    plot(t/60, qNormTrue - 1, 'r--', 'LineWidth', 1.2, 'DisplayName', 'True');
    
    ylabel('||q|| - 1', 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    title('Quaternion Normalization Check', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);
    
    set(gca, 'FontSize', 9);
    xlim([t(1)/60, t(end)/60]);
    ylim([-1e-10, 1e-10]);
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'MEKF_fig2_attitude_norm');
    end
    
    %% ------------------------------------------------------------------------
    % 3. GYRO BIAS ESTIMATION ERROR
    % -------------------------------------------------------------------------
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'MEKF - Gyro Bias Error', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [200, 200, 1370, 890]);
    
    axisLabels = {'X-axis [deg/h]', 'Y-axis [deg/h]', 'Z-axis [deg/h]'};
    
    for i = 1:3
        subplot(3,1,i);
        hold on; grid on;
        
        plot(t/60, biasErrorDegph(i,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Error');
        plot(t/60,  3*sigmaBias(i,:), 'r--', 'LineWidth', 1.2, 'DisplayName', '±3σ');
        plot(t/60, -3*sigmaBias(i,:), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        
        ylabel(axisLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 9);
        
        if i == 1
            title('Gyro Bias Estimation Error', 'FontSize', 13, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
        end
        
        set(gca, 'FontSize', 9);
        xlim([t(1)/60, t(end)/60]);
        
        rmsError = rms(biasErrorDegph(i,:));
        text(0.02, 0.95, sprintf('RMS = %.3f deg/h', rmsError), ...
             'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'white');
    end
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'MEKF_fig3_gyro_bias_error');
    end
    
    %% ------------------------------------------------------------------------
    % 4. ATTITUDE AND BIAS UNCERTAINTY EVOLUTION (3σ)
    % -------------------------------------------------------------------------
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'MEKF - Uncertainty Evolution', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [250, 250, 1370, 600]);
    
    subplot(2,1,1);
    hold on; grid on;
    
    plot(t/60, 3*sigmaAtt(1,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Roll');
    plot(t/60, 3*sigmaAtt(2,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Pitch');
    plot(t/60, 3*sigmaAtt(3,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Yaw');
    
    ylabel('3σ Uncertainty [deg]', 'FontSize', 10, 'FontWeight', 'bold');
    title('Attitude Uncertainty Evolution', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 10);
    
    set(gca, 'FontSize', 9);
    xlim([t(1)/60, t(end)/60]);
    
    subplot(2,1,2);
    hold on; grid on;
    
    plot(t/60, 3*sigmaBias(1,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'X-axis');
    plot(t/60, 3*sigmaBias(2,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Y-axis');
    plot(t/60, 3*sigmaBias(3,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Z-axis');
    
    ylabel('3σ Uncertainty [deg/h]', 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    title('Gyro Bias Uncertainty Evolution', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 10);
    
    set(gca, 'FontSize', 9);
    xlim([t(1)/60, t(end)/60]);
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'MEKF_fig4_uncertainty_evolution');
    end
    
    %% ------------------------------------------------------------------------
    % 5. NEES CONSISTENCY CHECK (CHI-SQUARED BOUNDS)
    % -------------------------------------------------------------------------
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'MEKF - NEES Consistency', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [300, 300, 1370, 600]);
    
    NEESatt = zeros(1, N);
    for k = 1:N
        % MEKF error state from quaternion error
        qw         = max(abs(qError(1,k)), 1e-10);
        deltaTheta = 2 * qError(2:4, k) / qw;
        
        % Add regularization to avoid numerical issues
        Preg  = P_hist(1:3, 1:3, k) + eye(3)*1e-12;
        Pinv  = inv(Preg);
        
        % NEES = δθ' * P^(-1) * δθ
        NEESatt(k) = deltaTheta' * Pinv * deltaTheta;
    end
    
    subplot(2,1,1);
    hold on; grid on;
    
    plot(t/60, NEESatt, 'b-', 'LineWidth', 1.5, 'DisplayName', 'NEES');
    
    % Chi-squared bounds (95% confidence interval for 3 DOF)
    r1 = chi2inv(0.025, 3);
    r2 = chi2inv(0.975, 3);
    
    plot(t/60, r2*ones(size(t)), 'r--', 'LineWidth', 1.2, 'DisplayName', '95% CI');
    plot(t/60, r1*ones(size(t)), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    
    ylabel('NEES', 'FontSize', 10, 'FontWeight', 'bold');
    title('Normalized Estimation Error Squared (NEES) - Attitude', ...
          'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);
    
    set(gca, 'FontSize', 9);
    xlim([t(1)/60, t(end)/60]);
    
    % Percentage inside bounds
    inside = sum(NEESatt > r1 & NEESatt < r2) / N * 100;
    text(0.02, 0.95, sprintf('%.1f%% inside 95%% CI', inside), ...
         'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'white', ...
         'VerticalAlignment', 'top');
    
    % Average NEES
    subplot(2,1,2);
    hold on; grid on;
    
    avgNEES = cumsum(NEESatt) ./ (1:N);
    plot(t/60, avgNEES, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Average NEES');
    plot(t/60, 3*ones(size(t)), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Expected (3 DOF)');
    
    ylabel('Average NEES', 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
    title('Cumulative Average NEES (should converge to DOF = 3)', ...
          'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);
    
    set(gca, 'FontSize', 9);
    xlim([t(1)/60, t(end)/60]);
    ylim([0 10]);
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'MEKF_fig5_NEES_consistency');
    end
    
    %% ------------------------------------------------------------------------
    % 6. QUATERNION COMPONENTS: TRUE VS ESTIMATED
    % -------------------------------------------------------------------------
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'MEKF - Quaternion Components', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [350, 350, 1370, 890]);
    
    labels = {'q_w (scalar)', 'q_x', 'q_y', 'q_z'};
    
    for i = 1:4
        subplot(4,1,i);
        hold on; grid on;
        
        plot(t/60, truth.qTrue(i,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'True');
        plot(t/60, qEst(i,:), 'b--', 'LineWidth', 1.2, 'DisplayName', 'Estimated');
        
        ylabel(labels{i}, 'FontSize', 10, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 9);
        
        if i == 1
            title('Quaternion Components: True vs Estimated', ...
                  'FontSize', 13, 'FontWeight', 'bold');
        end
        if i == 4
            xlabel('Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
        end
        
        set(gca, 'FontSize', 9);
        xlim([t(1)/60, t(end)/60]);
    end
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'MEKF_fig6_quaternion_components');
    end
    
    %% ------------------------------------------------------------------------
    % 7. ANGULAR VELOCITY: TRUE VS MEASURED
    % -------------------------------------------------------------------------
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'MEKF - Angular Velocity', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [400, 400, 1370, 800]);
    
    axisLabels = {'ω_x [deg/s]', 'ω_y [deg/s]', 'ω_z [deg/s]'};
    
    for i = 1:3
        subplot(3,1,i);
        hold on; grid on;
        
        plot(t/60, rad2deg(truth.omegaTrue(i,:)), 'r-', 'LineWidth', 1.5, ...
             'DisplayName', 'True');
        plot(t/60, rad2deg(imuMeas.gyro.omegaBody(i,:)), 'b--', 'MarkerSize', 1, ...
             'DisplayName', 'Measured');
        
        ylabel(axisLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 9);
        
        if i == 1
            title('Angular Velocity: True vs Measured (IMU)', ...
                  'FontSize', 13, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [min]', 'FontSize', 10, 'FontWeight', 'bold');
        end
        
        set(gca, 'FontSize', 9);
        xlim([t(1)/60, t(end)/60]);
    end
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'MEKF_fig7_angular_velocity');
    end
    
    %% ------------------------------------------------------------------------
    % 8. ERROR STATISTICS SUMMARY (BAR CHARTS)
    % -------------------------------------------------------------------------
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'MEKF - Error Statistics', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [450, 450, 1100, 600]);
    
    % Compute statistics
    stats.attMean  = mean(abs(attErrorDeg), 2);
    stats.attStd   = std(attErrorDeg, 0, 2);
    stats.attMax   = max(abs(attErrorDeg), [], 2);
    stats.attRms   = sqrt(mean(attErrorDeg.^2, 2));
    
    stats.biasMean = mean(abs(biasErrorDegph), 2);
    stats.biasStd  = std(biasErrorDegph, 0, 2);
    stats.biasMax  = max(abs(biasErrorDegph), [], 2);
    stats.biasRms  = sqrt(mean(biasErrorDegph.^2, 2));
    
    % Plot bar charts
    subplot(1,2,1);
    x     = 1:3;
    width = 0.2;
    
    bar(x - 1.5*width, stats.attMean, width, 'DisplayName', 'Mean');
    hold on;
    bar(x - 0.5*width, stats.attStd, width, 'DisplayName', 'Std Dev');
    bar(x + 0.5*width, stats.attRms, width, 'DisplayName', 'RMS');
    bar(x + 1.5*width, stats.attMax, width, 'DisplayName', 'Max');
    
    set(gca, 'XTick', 1:3, 'XTickLabel', {'Roll', 'Pitch', 'Yaw'}, 'FontSize', 10);
    ylabel('Error [deg]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Attitude Error Statistics', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);
    grid on;
    
    subplot(1,2,2);
    bar(x - 1.5*width, stats.biasMean, width, 'DisplayName', 'Mean');
    hold on;
    bar(x - 0.5*width, stats.biasStd, width, 'DisplayName', 'Std Dev');
    bar(x + 0.5*width, stats.biasRms, width, 'DisplayName', 'RMS');
    bar(x + 1.5*width, stats.biasMax, width, 'DisplayName', 'Max');
    
    set(gca, 'XTick', 1:3, 'XTickLabel', {'X', 'Y', 'Z'}, 'FontSize', 10);
    ylabel('Error [deg/h]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Gyro Bias Error Statistics', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);
    grid on;
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'MEKF_fig8_error_statistics');
    end
   
end
