function fig_handles = plotMAGResults(t, B_true_B, mag_meas, MAG, saveFlag)
%==========================================================================
% plotMAGResults - Generate comprehensive visualization of MAG simulation,
%                  calibration artifacts, and navigation/attitude deltas.
%
% INPUTS:
%   t             - Time vector [s], 1xN or Nx1
%   B_true_B      - True magnetic field in BODY frame [nT], 3xN
%   mag_meas      - MAG measurement struct:
%                     .B_meas     [nT], 3xN (measured, typically Sensor frame)
%                     .B_clean    [nT], 3xN (pre-quantization, Sensor frame)
%                     .B_true_S   [nT], 3xN (truth in Sensor frame)
%                     .bias_dyn   [nT], 3xN (dynamic bias in Sensor frame)
%                     .bias_total [nT], 3xN (total bias in Sensor frame)
%                     .B_det      [nT], 3xN (deterministic-only in Sensor frame)
%   MAG           - MAG parameter struct (optional fields used):
%                     .DCM_B2S_true    (Body to Sensor true mounting DCM)
%                     .resolution      [nT]
%                     .range           [nT]
%   saveFlag      - (Optional) If true, saves figures to Figures/MAG/
%                   Default: false
%
% OUTPUTS:
%   fig_handles   - Array of figure handles
%
% FIGURES GENERATED (as available):
%   1) Magnetic field error (components + norm) in BODY frame
%   2) Magnetic field direction error (angle) and magnitude error
%   3) Dynamic bias evolution (and total bias if available)
%   4) 3D B-vector cloud (measured vs truth) to visualize ellipsoid distortion
%      (useful for hard/soft-iron calibration validation)
%==========================================================================

    if nargin < 5
        saveFlag = false;
    end
    
    % Ensure time is row vector
    t = t(:).';
    N = numel(t);
    
    fig_handles = gobjects(4,1);
    Nfig = 0;
    
    % Prepare output directory if saving
    if saveFlag
        saveDir = fullfile(pwd, 'Figures', 'MAG');
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
            fprintf('\n✓ Created directory: %s\n', saveDir);
        end
        fprintf('\n--- Saving MAG figures ---\n');
    end
    
    % --- Frame harmonization ---
    % Truth in body frame is given: B_true_B (3xN)
    if ~isequal(size(B_true_B), [3, N])
        error('B_true_B must be 3xN with N = length(t).');
    end
    
    % Determine measured B and convert to body frame.
    B_meas_plot = MAG.DCM_B2S_true.' * mag_meas.B_meas;
    B_true_plot = B_true_B;
    
    %% ====================================================================
    % FIGURE 1: Magnetic Field Vector Error (Components + Norm)
    % =====================================================================
    if ~isempty(B_meas_plot) && ~isempty(B_true_plot)
        B_err = B_meas_plot - B_true_plot;
        B_err_norm = vecnorm(B_err, 2, 1);
        
        Nfig = Nfig + 1;
        fig_handles(Nfig) = figure('Name', 'MAG - B Vector Error', ...
            'Color', 'w', 'NumberTitle', 'off');
        
        axlbl = {'x','y','z'};
        for i = 1:3
            subplot(4,1,i);
            plot(t, B_err(i,:), 'b-', 'LineWidth', 1.2);
            grid on;
            ylabel(sprintf('\\DeltaB_%s [nT]', axlbl{i}), 'FontSize', 11, 'FontWeight', 'bold');
            if i == 1
                title(sprintf('Magnetic Field Error'), 'FontSize', 13, 'FontWeight', 'bold');
            end
            set(gca, 'FontSize', 10);
            xlim([t(1), t(end)]);
        end
        
        subplot(4,1,4);
        plot(t, B_err_norm, 'k-', 'LineWidth', 1.2);
        grid on;
        xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('\Delta||B|| [nT]', 'FontSize', 11, 'FontWeight', 'bold');
        set(gca, 'FontSize', 10);
        xlim([t(1), t(end)]);
        
        if saveFlag
            saveFigure(fig_handles(Nfig), saveDir, 'MAG_fig1_B_error');
        end
    end
    
    %% ====================================================================
    % FIGURE 2: Direction and Magnitude Error
    % =====================================================================
    u_true = B_true_plot ./ vecnorm(B_true_plot, 2, 1);
    u_meas = B_meas_plot ./ vecnorm(B_meas_plot, 2, 1);
    
    cang = sum(u_true .* u_meas, 1);
    cang = max(min(cang, 1), -1);
    ang_err = rad2deg(acos(cang));
    
    B_mag_true = vecnorm(B_true_plot, 2, 1);
    B_mag_meas = vecnorm(B_meas_plot, 2, 1);
    mag_err    = B_mag_meas - B_mag_true;
    
    Nfig = Nfig + 1;
    fig_handles(Nfig) = figure('Name', 'MAG - Direction & Magnitude Error', ...
        'Color', 'w', 'NumberTitle', 'off');
    
    subplot(2,1,1);
    plot(t, ang_err, 'm-', 'LineWidth', 1.2);
    grid on;
    ylabel('Angle [deg]', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('Field Direction Error'), 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontSize', 10);
    xlim([t(1), t(end)]);
    
    subplot(2,1,2);
    plot(t, mag_err, 'b-', 'LineWidth', 1.2);
    grid on;
    xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('||B|| [nT]', 'FontSize', 11, 'FontWeight', 'bold');
    set(gca, 'FontSize', 10);
    xlim([t(1), t(end)]);
    
    if saveFlag
        saveFigure(fig_handles(Nfig), saveDir, 'MAG_fig2_direction_magnitude_error');
    end
    
    %% ====================================================================
    % FIGURE 3: Dynamic Bias Evolution
    % =====================================================================
    Nfig = Nfig + 1;
    fig_handles(Nfig) = figure('Name', 'MAG - Bias Drift', ...
                               'Color', 'w', 'NumberTitle', 'off');

    bias_plot = MAG.DCM_B2S_true.' * mag_meas.bias_dyn;

    hold on;
    plot(t, bias_plot(1,:), 'r', 'LineWidth', 1.2, 'DisplayName', 'x');
    plot(t, bias_plot(2,:), 'g', 'LineWidth', 1.2, 'DisplayName', 'y');
    plot(t, bias_plot(3,:), 'b', 'LineWidth', 1.2, 'DisplayName', 'z');
    hold off;
    
    grid on;
    xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Dynamic Bias [nT]', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('MAG Dynamic Bias Random Walk'), 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northoutside', 'Orientation', 'horizontal', ...
           'FontSize', 10);
    set(gca, 'FontSize', 10);
    xlim([t(1), t(end)]);
    
    if saveFlag
        saveFigure(fig_handles(Nfig), saveDir, 'MAG_fig3_bias_drift');
    end
    
    %% ====================================================================
    % FIGURE 4: 3D B-Vector Cloud (Measured vs Truth) (Calibration Insight)
    % =====================================================================
    Nfig = Nfig + 1;
    fig_handles(Nfig) = figure('Name', 'MAG - 3D B Cloud', ...
        'Color', 'w', 'NumberTitle', 'off');
    
    % Downsample for plotting performance if needed
    idx = 1:N;
    if N > 8000
        idx = round(linspace(1, N, 8000));
    end
   
    hold on;
    plot3(B_true_plot(1,idx), B_true_plot(2,idx), B_true_plot(3,idx), ...
          'b-', 'DisplayName', 'Truth');
    plot3(B_meas_plot(1,idx), B_meas_plot(2,idx), B_meas_plot(3,idx), ...
          'r-', 'DisplayName', 'Measured');
    hold off;
    
    grid on; 
    axis equal;
    xlabel('B_x [nT]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('B_y [nT]', 'FontSize', 11, 'FontWeight', 'bold');
    zlabel('B_z [nT]', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('3D Magnetic Field Cloud'), 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northoutside', 'Orientation', 'horizontal', ...
           'FontSize', 10);
    set(gca, 'FontSize', 10);
    view(35, 25);
    
    if saveFlag
        saveFigure(fig_handles(Nfig), saveDir, 'MAG_fig4_3D_cloud');
    end
    
end