function fig_handles = plotIMUResults(t, omega_true_b, f_true_b, imu_meas, saveFlag)
%==========================================================================
% plotIMUResults: Generate quicklook figures for IMU truth vs measurements.
%
% Inputs:
%   t             - Time vector [s], 1xN
%   omega_true_b  - True angular rate in BODY frame [rad/s], 3xN
%   f_true_b      - True specific force in BODY frame [m/s^2], 3xN
%   imu_meas      - Struct with IMU measurements:
%                   .gyro_meas_b    [rad/s], 3xN
%                   .gyro_bias_dyn  [rad/s], 3xN
%                   .accel_meas_b   [m/s^2], 3xN
%                   .accel_bias_dyn [m/s^2], 3xN
%   saveFlag      - (Optional) Boolean. If true, saves figures to 'Figures/IMU'.
%                   Default: false
%
% Outputs:
%   fig_handles   - Array of figure handles
%
% Plots generated:
%   1) Gyro rate error and integrated angle error (approx orientation error)
%   2) Gyro dynamic bias drift
%   3) Accelerometer error and integrated velocity error
%   4) Accelerometer dynamic bias drift
%   5) Gyro Allan deviation
%   6) Accelerometer Allan deviation
%==========================================================================
    % Handle optional saveFlag
    if nargin < 5
        saveFlag = false;
    end
    
    Nfig = 0;
    fig_handles = gobjects(6,1);
    dt = t(2) - t(1);
    
    % Prepare output directory if saving
    if saveFlag
        saveDir = fullfile(pwd, 'Figures', 'IMU');
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
            fprintf('\n✓ Created directory: %s\n', saveDir);
        end
        fprintf('\n--- Saving IMU figures ---\n');
    end
    
    %% ------------------------------------------------------------------------
    % 1. GYRO RATE ERROR & INTEGRATED ANGLE ERROR (approx orientation error)
    % -------------------------------------------------------------------------
    gyro_err = imu_meas.gyro_meas_b - omega_true_b; % [rad/s], 3xN
    theta_err = cumsum(gyro_err, 2) * dt;           % [rad]
    Nfig = Nfig + 1;
    fig_handles(Nfig) = figure('Name', 'IMU - Gyro Errors', ...
                               'Color', 'w', 'NumberTitle', 'off', ...
                               'Position', [200, 100, 1370, 890]);
    axes_labels_rate  = {'\omega_x [deg/s]', '\omega_y [deg/s]', '\omega_z [deg/s]'};
    axes_labels_angle = {'\Delta\theta_x [deg]', '\Delta\theta_y [deg]', '\Delta\theta_z [deg]'};
    for i = 1:3
        subplot(3,2,2*i-1); % left: rate error
        plot(t, rad2deg(gyro_err(i,:)), 'r-', 'LineWidth', 1.0);
        grid on;
        ylabel(axes_labels_rate{i}, 'FontSize', 10, 'FontWeight', 'bold');
        if i == 1
            title('Gyro Rate Error', 'FontSize', 12, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [s]', 'FontSize', 10, 'FontWeight', 'bold');
        end
        set(gca, 'FontSize', 9);
        xlim([min(t(1)), max(t(end))]);
        
        subplot(3,2,2*i);   % right: integrated angle error
        plot(t, rad2deg(theta_err(i,:)), 'b-', 'LineWidth', 1.0);
        grid on;
        ylabel(axes_labels_angle{i}, 'FontSize', 10, 'FontWeight', 'bold');
        if i == 1
            title('Integrated Angular Error (small-angle approx.)', ...
                  'FontSize', 12, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [s]', 'FontSize', 10, 'FontWeight', 'bold');
        end
        set(gca, 'FontSize', 9);
        xlim([min(t(1)), max(t(end))]);
    end
    
    if saveFlag
        saveFigure(fig_handles(Nfig), saveDir, 'IMU_fig1_gyro_errors');
    end
    
    %% ------------------------------------------------------------------------
    % 2. GYRO DYNAMIC BIAS EVOLUTION
    % -------------------------------------------------------------------------
    Nfig = Nfig + 1;
    fig_handles(Nfig) = figure('Name', 'IMU - Gyro Dynamic Bias', ...
                               'Color', 'w', 'NumberTitle', 'off');
    plot(t, rad2deg(imu_meas.gyro_bias_dyn(1,:))*3600, 'r', 'LineWidth', 1.2, ...
        'DisplayName', 'x');
    hold on;
    plot(t, rad2deg(imu_meas.gyro_bias_dyn(2,:))*3600, 'g', 'LineWidth', 1.2, ...
        'DisplayName', 'y');
    plot(t, rad2deg(imu_meas.gyro_bias_dyn(3,:))*3600, 'b', 'LineWidth', 1.2, ...
        'DisplayName', 'z');
    grid on;
    xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Dynamic Bias [deg/h]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Gyro Dynamic Bias ARW', ...
          'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northoutside', 'Orientation', 'horizontal', ...
           'FontSize', 10);
    set(gca, 'FontSize', 10);
    xlim([min(t(1)), max(t(end))]);
    
    if saveFlag
        saveFigure(fig_handles(Nfig), saveDir, 'IMU_fig2_gyro_bias_drift');
    end
    
    %% ------------------------------------------------------------------------
    % 3. ACCELEROMETER ERROR & INTEGRATED VELOCITY ERROR
    % -------------------------------------------------------------------------
    accel_err = imu_meas.accel_meas_b - f_true_b;  % [m/s^2], 3xN
    vel_err = cumsum(accel_err, 2) * dt;           % [m/s]
    Nfig = Nfig + 1;
    fig_handles(Nfig) = figure('Name', 'IMU - Accelerometer Errors', ...
                               'Color', 'w', 'NumberTitle', 'off', ...
                               'Position', [200, 100, 1370, 890]);
    axes_labels_accel = {'a_x [m/s^2]', 'a_y [m/s^2]', 'a_z [m/s^2]'};
    axes_labels_vel   = {'\Deltav_x [m/s]', '\Deltav_y [m/s]', '\Deltav_z [m/s]'};
    for i = 1:3
        subplot(3,2,2*i-1); % left: acceleration error
        plot(t, accel_err(i,:), 'r-', 'LineWidth', 1.0);
        grid on;
        ylabel(axes_labels_accel{i}, 'FontSize', 10, 'FontWeight', 'bold');
        if i == 1
            title('Accelerometer Error', 'FontSize', 12, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [s]', 'FontSize', 10, 'FontWeight', 'bold');
        end
        set(gca, 'FontSize', 9);
        xlim([min(t(1)), max(t(end))]);
        
        subplot(3,2,2*i);   % right: integrated velocity error
        plot(t, vel_err(i,:), 'b-', 'LineWidth', 1.0);
        grid on;
        ylabel(axes_labels_vel{i}, 'FontSize', 10, 'FontWeight', 'bold');
        if i == 1
            title('Integrated Velocity Error', ...
                  'FontSize', 12, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [s]', 'FontSize', 10, 'FontWeight', 'bold');
        end
        set(gca, 'FontSize', 9);
        xlim([min(t(1)), max(t(end))]);
    end
    
    if saveFlag
        saveFigure(fig_handles(Nfig), saveDir, 'IMU_fig3_accel_errors');
    end
    
    %% ------------------------------------------------------------------------
    % 4. ACCELEROMETER DYNAMIC BIAS EVOLUTION
    % -------------------------------------------------------------------------
    Nfig = Nfig + 1;
    fig_handles(Nfig) = figure('Name', 'IMU - Accelerometer Dynamic Bias', ...
                               'Color', 'w', 'NumberTitle', 'off');
    plot(t, imu_meas.accel_bias_dyn(1,:), 'r', 'LineWidth', 1.2, ...
        'DisplayName', 'x');
    hold on;
    plot(t, imu_meas.accel_bias_dyn(2,:), 'g', 'LineWidth', 1.2, ...
        'DisplayName', 'y');
    plot(t, imu_meas.accel_bias_dyn(3,:), 'b', 'LineWidth', 1.2, ...
        'DisplayName', 'z');
    grid on;
    xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Dynamic Bias [m/s^2]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Accelerometer Dynamic Bias VRW', ...
          'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northoutside', 'Orientation', 'horizontal', ...
           'FontSize', 10);
    set(gca, 'FontSize', 10);
    xlim([min(t(1)), max(t(end))]);
    
    if saveFlag
        saveFigure(fig_handles(Nfig), saveDir, 'IMU_fig4_accel_bias_drift');
    end
    
    %% ------------------------------------------------------------------------
    % 5. ALLAN DEVIATION - GYRO (ORIENTATION ERROR PROXY)
    % -------------------------------------------------------------------------
    % Compute Allan Deviation for all three axes
    [taus_g_x, sigma_allan_g_x] = computeAllanDeviation(gyro_err(1,:), dt);
    [taus_g_y, sigma_allan_g_y] = computeAllanDeviation(gyro_err(2,:), dt);
    [taus_g_z, sigma_allan_g_z] = computeAllanDeviation(gyro_err(3,:), dt);
    
    Nfig = Nfig + 1;
    fig_handles(Nfig) = figure('Name', 'IMU - Gyro Allan Deviation', ...
                               'Color', 'w', 'NumberTitle', 'off');
    
    % Plot all three axes
    loglog(taus_g_x, rad2deg(sigma_allan_g_x), 'r-o', 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', 'r', 'DisplayName', 'x');
    hold on;
    loglog(taus_g_y, rad2deg(sigma_allan_g_y), 'g-o', 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', 'g', 'DisplayName', 'y');
    loglog(taus_g_z, rad2deg(sigma_allan_g_z), 'b-o', 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', 'b', 'DisplayName', 'z');
    hold off;
    
    grid on;
    xlabel('\tau [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('\sigma(\tau) [deg/s]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Gyro Allan Deviation', ...
         'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northoutside', 'Orientation', 'horizontal', ...
           'FontSize', 10);
    set(gca, 'FontSize', 10);
    xlim([min(taus_g_x(1)), max(taus_g_x(end))]);
    
    if saveFlag
        saveFigure(fig_handles(Nfig), saveDir, 'IMU_fig5_gyro_allan_dev');
    end
    
    %% ------------------------------------------------------------------------
    % 6. ALLAN DEVIATION - ACCELEROMETER
    % -------------------------------------------------------------------------
    % Compute Allan Deviation for all three axes
    [taus_a_x, sigma_allan_a_x] = computeAllanDeviation(accel_err(1,:), dt);
    [taus_a_y, sigma_allan_a_y] = computeAllanDeviation(accel_err(2,:), dt);
    [taus_a_z, sigma_allan_a_z] = computeAllanDeviation(accel_err(3,:), dt);
    
    Nfig = Nfig + 1;
    fig_handles(Nfig) = figure('Name', 'IMU - Accelerometer Allan Deviation', ...
                               'Color', 'w', 'NumberTitle', 'off');
    
    % Plot all three axes
    loglog(taus_a_x, sigma_allan_a_x, 'r-o', 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', 'r', 'DisplayName', 'x');
    hold on;
    loglog(taus_a_y, sigma_allan_a_y, 'g-o', 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', 'g', 'DisplayName', 'y');
    loglog(taus_a_z, sigma_allan_a_z, 'b-o', 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', 'b', 'DisplayName', 'z');
    hold off;
    
    grid on;
    xlabel('\tau [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('\sigma(\tau) [m/s^2]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Accelerometer Allan Deviation', ...
        'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northoutside', 'Orientation', 'horizontal', ...
           'FontSize', 10);
    set(gca, 'FontSize', 10);
    xlim([min(taus_a_x(1)), max(taus_a_x(end))]);
    
    if saveFlag
        saveFigure(fig_handles(Nfig), saveDir, 'IMU_fig6_accel_allan_dev');
    end
    
    % Trim unused handles
    fig_handles = fig_handles(1:Nfig);
    
end