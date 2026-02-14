function fig = plotIMUResults(t, omegaTrue_body, forceTrue_body, meas, saveFlag)
%==========================================================================
% plotIMUResults: Generate quicklook figures for IMU truth vs measurements.
%
% Inputs:
%   t              - Time vector                           [s], 1xN
%   omegaTrue_body - True angular rate in body frame   [rad/s], 3xN
%   forceTrue_body - True specific force in body frame [m/s^2], 3xN
%   imu_meas       - Struct with IMU measurements:
%                    .gyro.omegaBody                   [rad/s], 3xN
%                    .gyro.biasDym                     [rad/s], 3xN
%                    .accel.omegaBody                  [m/s^2], 3xN
%                    .accel.biasDyn                    [m/s^2], 3xN
%   saveFlag       - (Optional) Boolean. If true, saves figures to 'Figures/IMU'.
%                    Default: false
%
% Outputs:
%   fig            - Array of figure handles
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
    
    nFig = 0;
    fig = gobjects(6,1);
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
    % 1. GYRO RATE ERROR & INTEGRATED ANGLE ERROR
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'IMU - Gyro Errors', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [200, 100, 1370, 890]);

    gyroErr  = meas.gyro.omegaBody - omegaTrue_body; % [rad/s], 3xN
    thetaErr = cumsum(gyroErr, 2) * dt;              % [rad]
    
    axesLabels_rate  = {    '\omega_x [deg/s]',     '\omega_y [deg/s]',     '\omega_z [deg/s]'};
    axesLabels_angle = {'\Delta\theta_x [deg]', '\Delta\theta_y [deg]', '\Delta\theta_z [deg]'};
    for i = 1:3
        subplot(3,2,2*i-1);
        plot(t, rad2deg(gyroErr(i,:)), 'r-', 'LineWidth', 1.0);
        ylabel(axesLabels_rate{i}, 'FontSize', 10, 'FontWeight', 'bold');

        if i == 1
            title('Gyroscope Rate Error', 'FontSize', 12, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [s]', 'FontSize', 10, 'FontWeight', 'bold');
        end

        set(gca, 'FontSize', 9);
        grid on;
        xlim([min(t(1)), max(t(end))]);
        

        subplot(3,2,2*i);
        plot(t, rad2deg(thetaErr(i,:)), 'b-', 'LineWidth', 1.0);
        ylabel(axesLabels_angle{i}, 'FontSize', 10, 'FontWeight', 'bold');

        if i == 1
            title('Integrated Angular Error (small-angle approx.)', ...
                  'FontSize', 12, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [s]', 'FontSize', 10, 'FontWeight', 'bold');
        end

        set(gca, 'FontSize', 9);
        grid on;
        xlim([min(t(1)), max(t(end))]);
    end
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'IMU_fig1_gyro_errors');
    end
    
    %% ------------------------------------------------------------------------
    % 2. GYRO DYNAMIC BIAS EVOLUTION
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'IMU - Gyro Dynamic Bias', ...
                       'Color', 'w', 'NumberTitle', 'off');

    hold on;
    plot(t, rad2deg(meas.gyro.biasDyn(1,:))*3600, 'r', 'LineWidth', 1.2, ...
        'DisplayName', 'x');
    plot(t, rad2deg(meas.gyro.biasDyn(2,:))*3600, 'g', 'LineWidth', 1.2, ...
        'DisplayName', 'y');
    plot(t, rad2deg(meas.gyro.biasDyn(3,:))*3600, 'b', 'LineWidth', 1.2, ...
        'DisplayName', 'z');
    hold off;

    xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Dynamic Bias [deg/h]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Gyro Dynamic Bias ARW', ...
          'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northoutside', 'Orientation', 'horizontal', ...
           'FontSize', 10);

    set(gca, 'FontSize', 10);
    grid on;
    xlim([min(t(1)), max(t(end))]);
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'IMU_fig2_gyro_bias_drift');
    end
    
    %% ------------------------------------------------------------------------
    % 3. ACCELEROMETER ERROR & INTEGRATED VELOCITY ERROR
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'IMU - Accelerometer Errors', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [200, 100, 1370, 890]);

    accelErr = meas.accel.forceBody - forceTrue_body; % [m/s^2], 3xN
    velErr   = cumsum(accelErr, 2) * dt;              % [m/s]
    
    axesLabels_accel = {    'a_x [m/s^2]',     'a_y [m/s^2]',     'a_z [m/s^2]'};
    axesLabels_vel   = {'\Deltav_x [m/s]', '\Deltav_y [m/s]', '\Deltav_z [m/s]'};

    for i = 1:3
        subplot(3,2,2*i-1);
        plot(t, accelErr(i,:), 'r-', 'LineWidth', 1.0);
        ylabel(axesLabels_accel{i}, 'FontSize', 10, 'FontWeight', 'bold');

        if i == 1
            title('Accelerometer Error', 'FontSize', 12, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [s]', 'FontSize', 10, 'FontWeight', 'bold');
        end

        set(gca, 'FontSize', 9);
        grid on;
        xlim([min(t(1)), max(t(end))]);
        

        subplot(3,2,2*i);
        plot(t, velErr(i,:), 'b-', 'LineWidth', 1.0);
        ylabel(axesLabels_vel{i}, 'FontSize', 10, 'FontWeight', 'bold');
        if i == 1
            title('Integrated Velocity Error', ...
                  'FontSize', 12, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [s]', 'FontSize', 10, 'FontWeight', 'bold');
        end
        set(gca, 'FontSize', 9);
        grid on;
        xlim([min(t(1)), max(t(end))]);
    end
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'IMU_fig3_accel_errors');
    end
    
    %% ------------------------------------------------------------------------
    % 4. ACCELEROMETER DYNAMIC BIAS EVOLUTION
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'IMU - Accelerometer Dynamic Bias', ...
                       'Color', 'w', 'NumberTitle', 'off');

    hold on;
    plot(t, meas.accel.biasDyn(1,:), 'r', 'LineWidth', 1.2, ...
        'DisplayName', 'x');
    plot(t, meas.accel.biasDyn(2,:), 'g', 'LineWidth', 1.2, ...
        'DisplayName', 'y');
    plot(t, meas.accel.biasDyn(3,:), 'b', 'LineWidth', 1.2, ...
        'DisplayName', 'z');
    hold off;
    
    xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Dynamic Bias [m/s^2]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Accelerometer Dynamic Bias VRW', ...
          'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northoutside', 'Orientation', 'horizontal', ...
           'FontSize', 10);

    set(gca, 'FontSize', 10);
    grid on;
    xlim([min(t(1)), max(t(end))]);
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'IMU_fig4_accel_bias_drift');
    end
    
    %% ------------------------------------------------------------------------
    % 5. ALLAN DEVIATION - GYRO (ORIENTATION ERROR PROXY)
    % -------------------------------------------------------------------------

    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'IMU - Gyro Allan Deviation', ...
                       'Color', 'w', 'NumberTitle', 'off');

    % Compute Allan Deviation for all three axes
    [gyroTau_x, gyroADEV_x] = computeAllanDeviation(gyroErr(1,:), dt);
    [gyroTau_y, gyroADEV_y] = computeAllanDeviation(gyroErr(2,:), dt);
    [gyroTau_z, gyroADEV_z] = computeAllanDeviation(gyroErr(3,:), dt);
    
    % Plot all three axes
    loglog(gyroTau_x, rad2deg(gyroADEV_x), 'r-o', 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', 'r', 'DisplayName', 'x');
    hold on;
    loglog(gyroTau_y, rad2deg(gyroADEV_y), 'g-o', 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', 'g', 'DisplayName', 'y');
    loglog(gyroTau_z, rad2deg(gyroADEV_z), 'b-o', 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', 'b', 'DisplayName', 'z');
    hold off;
    
    xlabel('\tau [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('\sigma(\tau) [deg/s]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Gyro Allan Deviation', ...
         'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northoutside', 'Orientation', 'horizontal', ...
           'FontSize', 10);

    set(gca, 'FontSize', 10);
    grid on;
    xlim([min(gyroTau_x(1)), max(gyroTau_x(end))]);
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'IMU_fig5_gyro_allan_dev');
    end
    
    %% ------------------------------------------------------------------------
    % 6. ALLAN DEVIATION - ACCELEROMETER
    % -------------------------------------------------------------------------
    
    nFig = nFig + 1;
    fig(nFig) = figure('Name', 'IMU - Accelerometer Allan Deviation', ...
                               'Color', 'w', 'NumberTitle', 'off');
    % Compute Allan Deviation for all three axes
    [accelTau_x, accelADEV_x] = computeAllanDeviation(accelErr(1,:), dt);
    [accelTau_y, accelADEV_y] = computeAllanDeviation(accelErr(2,:), dt);
    [accelTau_z, accelADEV_z] = computeAllanDeviation(accelErr(3,:), dt);
    
    % Plot all three axes
    loglog(accelTau_x, accelADEV_x, 'r-o', 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', 'r', 'DisplayName', 'x');
    hold on;
    loglog(accelTau_y, accelADEV_y, 'g-o', 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', 'g', 'DisplayName', 'y');
    loglog(accelTau_z, accelADEV_z, 'b-o', 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', 'b', 'DisplayName', 'z');
    hold off;
    
    xlabel('\tau [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('\sigma(\tau) [m/s^2]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Accelerometer Allan Deviation', ...
        'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northoutside', 'Orientation', 'horizontal', ...
           'FontSize', 10);

    set(gca, 'FontSize', 10);
    grid on;
    xlim([min(accelTau_x(1)), max(accelTau_x(end))]);
    
    if saveFlag
        saveFigure(fig(nFig), saveDir, 'IMU_fig6_accel_allan_dev');
    end

end