function plotEKFResults(t, truth, qEst, biasEst, P_hist, imuMeas, saveFlag)
%==========================================================================
% plotEKFResults - Comprehensive MEKF performance analysis and visualization
%
% INPUTS:
%   t        - Time vector [s], 1xN
%   truth    - Truth structure with fields:
%              .qTrue      - True quaternion [4xN]
%              .omegaTrue  - True angular velocity [rad/s, 3xN]
%              .rECI, .vECI, .B_ECI, etc.
%   qEst     - Estimated quaternion [4xN]
%   biasEst  - Estimated gyro bias [rad/s, 3xN]
%   P_hist   - Covariance history [6x6xN]
%   imuMeas  - IMU measurements structure
%   saveFlag - (OPTIONAL) Boolean flag to save figures [default: false]
%              If true, saves all figures in './results/MEKF/' directory
%
% OUTPUTS:
%   Multiple figures with performance analysis (displayed and optionally saved)
%==========================================================================

    fprintf('\n=== Generating MEKF Performance Plots ===\n');
    
    %% Handle optional saveFlag argument
    if nargin < 7
        saveFlag = false;  % Default: don't save figures
    end
    
    %% Create output directory if saving is enabled
    if saveFlag
        figDir = fullfile(pwd, 'Figures', 'MEKF');
        if ~exist(figDir, 'dir')
            mkdir(figDir);
            fprintf('Created directory: %s\n', figDir);
        end
    end
    
    N = length(t);
    
    %% ===================================================================
    % COMPUTE ERRORS
    % ====================================================================
    
    % 1. Attitude Error (Quaternion Error)
    qError = zeros(4, N);
    attError_deg = zeros(3, N);  % Euler angle representation [yaw; pitch; roll] in [deg]
    
    for k = 1:N
        % Quaternion error: q_err = q_true^{-1} ⊗ q_est
        qError(:,k) = quatmultiply(quatinv(truth.qTrue(:,k)), qEst(:,k));
        
        % Convert quaternion error to DCM, then extract Euler angles (ZYX)
        DCM_error = quat2dcm(qError(:,k));
        euler_error_rad = dcm2euler_ZYX(DCM_error);  % [yaw; pitch; roll]
        attError_deg(:,k) = rad2deg(euler_error_rad);
    end
    
    % 2. Gyro Bias Error (assume true bias is zero for simulation)
    biasTrue = zeros(3, N);  % If you have true bias, replace this
    biasError_degph = rad2deg(biasEst - biasTrue) * 3600;  % [deg/h]
    
    % 3. Extract Uncertainty (1σ) from Covariance
    sigma_att = zeros(3, N);    % Attitude uncertainty [deg]
    sigma_bias = zeros(3, N);   % Bias uncertainty [deg/h]
    
    for k = 1:N
        sigma_att(:,k)  = rad2deg(sqrt(diag(P_hist(1:3, 1:3, k))));
        sigma_bias(:,k) = rad2deg(sqrt(diag(P_hist(4:6, 4:6, k)))) * 3600;
    end
    
    %% ===================================================================
    % FIGURE 1: ATTITUDE ERROR (3-AXIS) WITH UNCERTAINTY BOUNDS
    % ====================================================================
    fig1 = figure('Name', 'MEKF - Attitude Error', 'Position', [100 100 1200 800]);
    
    axisLabels = {'Roll (X)', 'Pitch (Y)', 'Yaw (Z)'};
    
    for i = 1:3
        subplot(3,1,i);
        hold on; grid on;
        
        % Plot error
        plot(t/60, attError_deg(i,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Error');
        
        % Plot ±3σ bounds
        plot(t/60, 3*sigma_att(i,:), 'r--', 'LineWidth', 1.2, 'DisplayName', '±3σ');
        plot(t/60, -3*sigma_att(i,:), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        
        ylabel([axisLabels{i}, ' [deg]'], 'FontSize', 10);
        legend('Location', 'best', 'FontSize', 9);
        
        if i == 1
            title('Attitude Estimation Error (Euler Angles)', 'FontSize', 12, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [min]', 'FontSize', 10);
        end
        
        % Compute RMS error
        rms_error = rms(attError_deg(i,:));
        text(0.02, 0.95, sprintf('RMS = %.3f°', rms_error), ...
             'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'white');
    end
    
    if saveFlag
        saveFigure(fig1, figDir, 'MEKF_fig1_attitude_error');
    end
    
    %% ===================================================================
    % FIGURE 2: ATTITUDE ERROR NORM & QUATERNION COMPONENTS
    % ====================================================================
    fig2 = figure('Name', 'MEKF - Attitude Error Norm', 'Position', [150 150 1200 600]);
    
    % Subplot 1: Attitude Error Norm
    subplot(2,1,1);
    hold on; grid on;
    
    % Total attitude error (angle of rotation from quaternion error)
    attErrorNorm = zeros(1, N);
    for k = 1:N
        attErrorNorm(k) = 2 * acos(min(abs(qError(1,k)), 1));  % [rad]
    end
    
    plot(t/60, rad2deg(attErrorNorm), 'b-', 'LineWidth', 1.5);
    ylabel('Attitude Error [deg]', 'FontSize', 10);
    title('Total Attitude Error Magnitude', 'FontSize', 12, 'FontWeight', 'bold');
    
    % Statistics
    mean_err = mean(rad2deg(attErrorNorm));
    max_err = max(rad2deg(attErrorNorm));
    text(0.02, 0.95, sprintf('Mean = %.3f°\nMax = %.3f°', mean_err, max_err), ...
         'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'white', ...
         'VerticalAlignment', 'top');
    
    % Subplot 2: Quaternion Normalization Check
    subplot(2,1,2);
    hold on; grid on;
    
    qNorm_est = vecnorm(qEst, 2, 1);
    qNorm_true = vecnorm(truth.qTrue, 2, 1);
    
    plot(t/60, qNorm_est - 1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Estimated');
    plot(t/60, qNorm_true - 1, 'r--', 'LineWidth', 1.2, 'DisplayName', 'True');
    
    ylabel('||q|| - 1', 'FontSize', 10);
    xlabel('Time [min]', 'FontSize', 10);
    title('Quaternion Normalization Check', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best');
    ylim([-1e-10, 1e-10]);
    
    if saveFlag
        saveFigure(fig2, figDir, 'MEKF_fig2_attitude_norm');
    end
    
    %% ===================================================================
    % FIGURE 3: GYRO BIAS ESTIMATION ERROR
    % ====================================================================
    fig3 = figure('Name', 'MEKF - Gyro Bias Error', 'Position', [200 200 1200 800]);
    
    axisLabels = {'X-axis', 'Y-axis', 'Z-axis'};
    
    for i = 1:3
        subplot(3,1,i);
        hold on; grid on;
        
        % Plot bias error
        plot(t/60, biasError_degph(i,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Error');
        
        % Plot ±3σ bounds
        plot(t/60, 3*sigma_bias(i,:), 'r--', 'LineWidth', 1.2, 'DisplayName', '±3σ');
        plot(t/60, -3*sigma_bias(i,:), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        
        ylabel([axisLabels{i}, ' [deg/h]'], 'FontSize', 10);
        legend('Location', 'best', 'FontSize', 9);
        
        if i == 1
            title('Gyro Bias Estimation Error', 'FontSize', 12, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [min]', 'FontSize', 10);
        end
        
        % Compute RMS error
        rms_error = rms(biasError_degph(i,:));
        text(0.02, 0.95, sprintf('RMS = %.3f deg/h', rms_error), ...
             'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'white');
    end
    
    if saveFlag
        saveFigure(fig3, figDir, 'MEKF_fig3_gyro_bias_error');
    end
    
    %% ===================================================================
    % FIGURE 4: ATTITUDE UNCERTAINTY EVOLUTION (3σ)
    % ====================================================================
    fig4 = figure('Name', 'MEKF - Uncertainty Evolution', 'Position', [250 250 1200 600]);
    
    % Subplot 1: Attitude Uncertainty
    subplot(2,1,1);
    hold on; grid on;
    
    plot(t/60, 3*sigma_att(1,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Roll');
    plot(t/60, 3*sigma_att(2,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Pitch');
    plot(t/60, 3*sigma_att(3,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Yaw');
    
    ylabel('3σ Uncertainty [deg]', 'FontSize', 10);
    title('Attitude Uncertainty Evolution', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best');
    
    % Subplot 2: Bias Uncertainty
    subplot(2,1,2);
    hold on; grid on;
    
    plot(t/60, 3*sigma_bias(1,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'X-axis');
    plot(t/60, 3*sigma_bias(2,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Y-axis');
    plot(t/60, 3*sigma_bias(3,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Z-axis');
    
    ylabel('3σ Uncertainty [deg/h]', 'FontSize', 10);
    xlabel('Time [min]', 'FontSize', 10);
    title('Gyro Bias Uncertainty Evolution', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best');
    
    if saveFlag
        saveFigure(fig4, figDir, 'MEKF_fig4_uncertainty_evolution');
    end
    
    %% ===================================================================
    % FIGURE 5: NEES (Normalized Estimation Error Squared) - CONSISTENCY CHECK
    % ====================================================================
    fig5 = figure('Name', 'MEKF - NEES Consistency', 'Position', [300 300 1200 600]);
    
    % Compute NEES for attitude
    NEES_att = zeros(1, N);
    for k = 1:N
        % Error state (attitude part only, small angle)
        delta_theta = 2 * qError(2:4, k);  % Small-angle approximation
        
        % NEES = e' * P^{-1} * e
        P_inv = inv(P_hist(1:3, 1:3, k));
        NEES_att(k) = delta_theta' * P_inv * delta_theta;
    end
    
    % Plot NEES
    subplot(2,1,1);
    hold on; grid on;
    
    plot(t/60, NEES_att, 'b-', 'LineWidth', 1.5, 'DisplayName', 'NEES');
    
    % Chi-squared bounds (95% confidence interval for 3 DOF)
    r1 = chi2inv(0.025, 3);
    r2 = chi2inv(0.975, 3);
    
    plot(t/60, r2*ones(size(t)), 'r--', 'LineWidth', 1.2, 'DisplayName', '95% CI');
    plot(t/60, r1*ones(size(t)), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    
    ylabel('NEES', 'FontSize', 10);
    title('Normalized Estimation Error Squared (NEES) - Attitude', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best');
    
    % Percentage inside bounds
    inside = sum(NEES_att > r1 & NEES_att < r2) / N * 100;
    text(0.02, 0.95, sprintf('%.1f%% inside 95%% CI', inside), ...
         'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white', ...
         'VerticalAlignment', 'top');
    
    % Average NEES
    subplot(2,1,2);
    hold on; grid on;
    
    avgNEES = cumsum(NEES_att) ./ (1:N);
    plot(t/60, avgNEES, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Average NEES');
    plot(t/60, 3*ones(size(t)), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Expected (3 DOF)');
    
    ylabel('Average NEES', 'FontSize', 10);
    xlabel('Time [min]', 'FontSize', 10);
    title('Cumulative Average NEES (should converge to DOF = 3)', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best');
    ylim([0 10]);
    
    if saveFlag
        saveFigure(fig5, figDir, 'MEKF_fig5_NEES_consistency');
    end
    
    %% ===================================================================
    % FIGURE 6: QUATERNION COMPONENTS COMPARISON
    % ====================================================================
    fig6 = figure('Name', 'MEKF - Quaternion Components', 'Position', [350 350 1200 800]);
    
    labels = {'q_w (scalar)', 'q_x', 'q_y', 'q_z'};
    
    for i = 1:4
        subplot(4,1,i);
        hold on; grid on;
        
        plot(t/60, truth.qTrue(i,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'True');
        plot(t/60, qEst(i,:), 'b--', 'LineWidth', 1.2, 'DisplayName', 'Estimated');
        
        ylabel(labels{i}, 'FontSize', 10);
        legend('Location', 'best', 'FontSize', 9);
        
        if i == 1
            title('Quaternion Components: True vs Estimated', 'FontSize', 12, 'FontWeight', 'bold');
        end
        if i == 4
            xlabel('Time [min]', 'FontSize', 10);
        end
    end
    
    if saveFlag
        saveFigure(fig6, figDir, 'MEKF_fig6_quaternion_components');
    end
    
    %% ===================================================================
    % FIGURE 7: ANGULAR VELOCITY COMPARISON
    % ====================================================================
    fig7 = figure('Name', 'MEKF - Angular Velocity', 'Position', [400 400 1200 700]);
    
    axisLabels = {'ω_x', 'ω_y', 'ω_z'};
    
    for i = 1:3
        subplot(3,1,i);
        hold on; grid on;
        
        plot(t/60, rad2deg(truth.omegaTrue(i,:)), 'r-', 'LineWidth', 1.5, 'DisplayName', 'True');
        plot(t/60, rad2deg(imuMeas.gyro.omegaBody(i,:)), 'b.', 'MarkerSize', 1, 'DisplayName', 'Measured (noisy)');
        
        ylabel([axisLabels{i}, ' [deg/s]'], 'FontSize', 10);
        legend('Location', 'best', 'FontSize', 9);
        
        if i == 1
            title('Angular Velocity: True vs Measured (IMU)', 'FontSize', 12, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [min]', 'FontSize', 10);
        end
    end
    
    if saveFlag
        saveFigure(fig7, figDir, 'MEKF_fig7_angular_velocity');
    end
    
    %% ===================================================================
    % FIGURE 8: ESTIMATION ERROR STATISTICS SUMMARY
    % ====================================================================
    fig8 = figure('Name', 'MEKF - Error Statistics', 'Position', [450 450 900 600]);
    
    % Compute statistics
    stats = struct();
    stats.att_mean = mean(abs(attError_deg), 2);
    stats.att_std = std(attError_deg, 0, 2);
    stats.att_max = max(abs(attError_deg), [], 2);
    stats.att_rms = sqrt(mean(attError_deg.^2, 2));
    
    stats.bias_mean = mean(abs(biasError_degph), 2);
    stats.bias_std = std(biasError_degph, 0, 2);
    stats.bias_max = max(abs(biasError_degph), [], 2);
    stats.bias_rms = sqrt(mean(biasError_degph.^2, 2));
    
    % Plot bar charts
    subplot(1,2,1);
    x = 1:3;
    width = 0.2;
    
    bar(x - 1.5*width, stats.att_mean, width, 'DisplayName', 'Mean');
    hold on;
    bar(x - 0.5*width, stats.att_std, width, 'DisplayName', 'Std Dev');
    bar(x + 0.5*width, stats.att_rms, width, 'DisplayName', 'RMS');
    bar(x + 1.5*width, stats.att_max, width, 'DisplayName', 'Max');
    
    set(gca, 'XTick', 1:3, 'XTickLabel', {'Roll', 'Pitch', 'Yaw'});
    ylabel('Error [deg]', 'FontSize', 10);
    title('Attitude Error Statistics', 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);
    grid on;
    
    subplot(1,2,2);
    bar(x - 1.5*width, stats.bias_mean, width, 'DisplayName', 'Mean');
    hold on;
    bar(x - 0.5*width, stats.bias_std, width, 'DisplayName', 'Std Dev');
    bar(x + 0.5*width, stats.bias_rms, width, 'DisplayName', 'RMS');
    bar(x + 1.5*width, stats.bias_max, width, 'DisplayName', 'Max');
    
    set(gca, 'XTick', 1:3, 'XTickLabel', {'X', 'Y', 'Z'});
    ylabel('Error [deg/h]', 'FontSize', 10);
    title('Gyro Bias Error Statistics', 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);
    grid on;
    
    if saveFlag
        saveFigure(fig8, figDir, 'MEKF_fig8_error_statistics');
    end
    
    %% ===================================================================
    % PRINT SUMMARY STATISTICS TO CONSOLE
    % ====================================================================
    fprintf('\n========================================\n');
    fprintf('MEKF PERFORMANCE SUMMARY\n');
    fprintf('========================================\n');
    
    fprintf('\nAttitude Estimation Error:\n');
    fprintf('  Roll  - RMS: %.4f°, Max: %.4f°\n', stats.att_rms(1), stats.att_max(1));
    fprintf('  Pitch - RMS: %.4f°, Max: %.4f°\n', stats.att_rms(2), stats.att_max(2));
    fprintf('  Yaw   - RMS: %.4f°, Max: %.4f°\n', stats.att_rms(3), stats.att_max(3));
    
    fprintf('\nGyro Bias Estimation Error:\n');
    fprintf('  X-axis - RMS: %.4f deg/h, Max: %.4f deg/h\n', stats.bias_rms(1), stats.bias_max(1));
    fprintf('  Y-axis - RMS: %.4f deg/h, Max: %.4f deg/h\n', stats.bias_rms(2), stats.bias_max(2));
    fprintf('  Z-axis - RMS: %.4f deg/h, Max: %.4f deg/h\n', stats.bias_rms(3), stats.bias_max(3));
    
    fprintf('\nFilter Consistency (NEES):\n');
    fprintf('  %.1f%% of samples inside 95%% confidence bounds\n', inside);
    fprintf('  Average NEES: %.2f (expected = 3 for 3 DOF)\n', mean(NEES_att));
    
    fprintf('\n========================================\n');
    
    if saveFlag
        fprintf('\n=== All figures saved to: %s ===\n', figDir);
    end
    
    fprintf('=== Plotting Complete ===\n\n');
end
