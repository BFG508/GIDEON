function fig_handles = plotIMUResults(t, omega_true_b, f_true_b, imu_meas, use_accel_truth)
%==========================================================================
% plotIMUResults: Generate quicklook figures for IMU truth vs measurements.
%
% Inputs:
%   t             - Time vector [s], 1xN
%   omega_true_b  - True angular rate in BODY frame [rad/s], 3xN
%   f_true_b      - True specific force in BODY frame [m/s^2], 3xN or []
%   imu_meas      - Struct with IMU measurements (from generateIMUMeasurements)
%                   .gyro_meas_b    [rad/s], 3xN
%                   .gyro_bias_dyn  [rad/s], 3xN
%                   .accel_meas_b   [m/s^2], 3xN (optional)
%                   .accel_bias_dyn [m/s^2], 3xN (optional)
%   use_accel_truth - Boolean flag indicating if accelerometer truth is used
%
% Outputs:
%   fig_handles   - Array of figure handles
%
% Plots generated:
%   1) Gyro true vs measured angular rate (3 subplots)
%   2) Gyro rate error and integrated angle error (approx orientation error)
%   3) Gyro dynamic bias drift
%   4) (Optional) Accel true vs measured specific force
%   5) (Optional) Accel error and bias drift
%   6) Gyro Allan deviation (orientation error proxy, 1 axis)
%   7) (Optional) Accelerometer Allan deviation (1 axis)
%==========================================================================

Nfig = 0;
fig_handles = gobjects(7,1); % maximum expected
dt = t(2) - t(1);

%% ------------------------------------------------------------------------
% 1. GYRO: TRUE vs MEASURED RATES
% -------------------------------------------------------------------------
Nfig = Nfig + 1;
fig_handles(Nfig) = figure('Name', 'IMU - Gyro Measurements', ...
                           'Color', 'w', 'NumberTitle', 'off');

axes_labels = {'\omega_x [deg/s]', '\omega_y [deg/s]', '\omega_z [deg/s]'};

for i = 1:3
    subplot(3,1,i);
    plot(t, rad2deg(omega_true_b(i,:)), 'k-', 'LineWidth', 1.5, ...
        'DisplayName', 'True');
    hold on;
    plot(t, rad2deg(imu_meas.gyro_meas_b(i,:)), 'b-', 'LineWidth', 1.0, ...
        'DisplayName', 'Measured');
    grid on;
    ylabel(axes_labels{i}, 'FontSize', 11, 'FontWeight', 'bold');
    if i == 1
        title('Gyro Angular Rate: True vs Measured', ...
            'FontSize', 13, 'FontWeight', 'bold');
    end
    if i == 3
        xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold');
    end
    legend('Location', 'best');
    set(gca, 'FontSize', 10);
end

%% ------------------------------------------------------------------------
% 2. GYRO RATE ERROR & INTEGRATED ANGLE ERROR (approx orientation error)
% -------------------------------------------------------------------------
gyro_err = imu_meas.gyro_meas_b - omega_true_b;    % [rad/s], 3xN
theta_err = cumsum(gyro_err, 2) * dt;              % [rad]

Nfig = Nfig + 1;
fig_handles(Nfig) = figure('Name', 'IMU - Gyro Errors', ...
                           'Color', 'w', 'NumberTitle', 'off');

axes_labels_rate  = {'Error \omega_x [deg/s]', 'Error \omega_y [deg/s]', 'Error \omega_z [deg/s]'};
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
end

%% ------------------------------------------------------------------------
% 3. GYRO DYNAMIC BIAS EVOLUTION
% -------------------------------------------------------------------------
Nfig = Nfig + 1;
fig_handles(Nfig) = figure('Name', 'IMU - Gyro Dynamic Bias', ...
                           'Color', 'w', 'NumberTitle', 'off');

plot(t, rad2deg(imu_meas.gyro_bias_dyn(1,:))*3600, 'r', 'LineWidth', 1.2, ...
    'DisplayName', 'Bias_x');
hold on;
plot(t, rad2deg(imu_meas.gyro_bias_dyn(2,:))*3600, 'g', 'LineWidth', 1.2, ...
    'DisplayName', 'Bias_y');
plot(t, rad2deg(imu_meas.gyro_bias_dyn(3,:))*3600, 'b', 'LineWidth', 1.2, ...
    'DisplayName', 'Bias_z');
grid on;
xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Dynamic Bias [deg/h]', 'FontSize', 11, 'FontWeight', 'bold');
title('Gyro Dynamic Bias Random Walk', ...
    'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best');
set(gca, 'FontSize', 10);

%% ------------------------------------------------------------------------
% 4. OPTIONAL: ACCELEROMETER QUICKLOOK (if enabled)
% -------------------------------------------------------------------------
if use_accel_truth && ~isempty(f_true_b) && isfield(imu_meas, 'accel_meas_b')
    % 4.1 True vs measured specific force
    Nfig = Nfig + 1;
    fig_handles(Nfig) = figure('Name', 'IMU - Accelerometer Measurements', ...
                               'Color', 'w', 'NumberTitle', 'off');

    axes_labels_acc = {'a_x [m/s^2]', 'a_y [m/s^2]', 'a_z [m/s^2]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, f_true_b(i,:), 'k-', 'LineWidth', 1.5, ...
            'DisplayName', 'True');
        hold on;
        plot(t, imu_meas.accel_meas_b(i,:), 'b-', 'LineWidth', 1.0, ...
            'DisplayName', 'Measured');
        grid on;
        ylabel(axes_labels_acc{i}, 'FontSize', 11, 'FontWeight', 'bold');
        if i == 1
            title('Accelerometer Specific Force: True vs Measured', ...
                'FontSize', 13, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold');
        end
        legend('Location', 'best');
        set(gca, 'FontSize', 10);
    end

    % 4.2 Accel error & bias drift
    accel_err = imu_meas.accel_meas_b - f_true_b;    % [m/s^2]
    Nfig = Nfig + 1;
    fig_handles(Nfig) = figure('Name', 'IMU - Accelerometer Errors', ...
                               'Color', 'w', 'NumberTitle', 'off');

    axes_labels_acc_err  = {'Error a_x [m/s^2]', 'Error a_y [m/s^2]', 'Error a_z [m/s^2]'};

    for i = 1:3
        subplot(3,2,2*i-1);
        plot(t, accel_err(i,:), 'r-', 'LineWidth', 1.0);
        grid on;
        ylabel(axes_labels_acc_err{i}, 'FontSize', 10, 'FontWeight', 'bold');
        if i == 1
            title('Accelerometer Error', ...
                  'FontSize', 12, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [s]', 'FontSize', 10, 'FontWeight', 'bold');
        end
        set(gca, 'FontSize', 9);

        subplot(3,2,2*i);
        plot(t, imu_meas.accel_bias_dyn(i,:), 'b-', 'LineWidth', 1.0);
        grid on;
        ylabel('Bias [m/s^2]', 'FontSize', 10, 'FontWeight', 'bold');
        if i == 1
            title('Accelerometer Dynamic Bias', ...
                  'FontSize', 12, 'FontWeight', 'bold');
        end
        if i == 3
            xlabel('Time [s]', 'FontSize', 10, 'FontWeight', 'bold');
        end
        set(gca, 'FontSize', 9);
    end
else
    accel_err = [];
end

%% ------------------------------------------------------------------------
% 5. ALLAN DEVIATION - GYRO (ORIENTATION ERROR PROXY)
% -------------------------------------------------------------------------
% We use rate error of X-axis as representative (constant truth recommended
% for rigorous analysis; here this is illustrative).

[taus_g, sigma_allan_g] = computeAllanDeviation(gyro_err(1,:), dt);

Nfig = Nfig + 1;
fig_handles(Nfig) = figure('Name', 'IMU - Gyro Allan Deviation (X-axis)', ...
                           'Color', 'w', 'NumberTitle', 'off');

loglog(taus_g, rad2deg(sigma_allan_g), 'b-o', 'LineWidth', 1.5, ...
    'MarkerSize', 5, 'MarkerFaceColor', 'b');
grid on;
xlabel('\tau [s]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('\sigma(\tau) [deg/s]', 'FontSize', 11, 'FontWeight', 'bold');
title('Gyro Allan Deviation (Rate Error, X-axis)', ...
    'FontSize', 13, 'FontWeight', 'bold');
set(gca, 'FontSize', 10);

%% ------------------------------------------------------------------------
% 6. ALLAN DEVIATION - ACCELEROMETER (if available)
% -------------------------------------------------------------------------
if ~isempty(accel_err)
    [taus_a, sigma_allan_a] = computeAllanDeviation(accel_err(1,:), dt);

    Nfig = Nfig + 1;
    fig_handles(Nfig) = figure('Name', 'IMU - Accelerometer Allan Deviation (X-axis)', ...
                               'Color', 'w', 'NumberTitle', 'off');

    loglog(taus_a, sigma_allan_a, 'r-o', 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', 'r');
    grid on;
    xlabel('\tau [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('\sigma(\tau) [m/s^2]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Accelerometer Allan Deviation (Error, X-axis)', ...
        'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontSize', 10);
end

% Trim unused handles
fig_handles = fig_handles(1:Nfig);

end

%==========================================================================
% HELPER: Compute Allan Deviation for a 1D time series
%==========================================================================
function [taus, sigma_allan] = computeAllanDeviation(x, dt)
% x: 1xN or Nx1 time series (zero-mean error recommended)
% dt: sample time [s]

x = x(:);              % column
N = length(x);

% Choose cluster sizes m (number of samples per averaging interval)
% Log-spaced from 1 sample to N/10 approximately
max_m = floor(N/10);
if max_m < 2
    taus = dt;
    sigma_allan = std(x);
    return;
end

m_list = unique(round(logspace(0, log10(max_m), 30)));
taus = m_list * dt;
sigma_allan = zeros(size(taus));

for k = 1:length(m_list)
    m = m_list(k);
    % Number of non-overlapping clusters of length m
    K = floor(N / m);
    if K < 2
        sigma_allan(k) = NaN;
        continue;
    end

    % Compute cluster averages
    y = zeros(K,1);
    for i = 1:K
        idx_start = (i-1)*m + 1;
        idx_end   = i*m;
        y(i) = mean(x(idx_start:idx_end));
    end

    % Allan variance
    dy = diff(y);
    sigma2 = 0.5 * mean(dy.^2);
    sigma_allan(k) = sqrt(sigma2);
end

% Remove NaNs (if any)
valid = ~isnan(sigma_allan);
taus = taus(valid);
sigma_allan = sigma_allan(valid);

end
