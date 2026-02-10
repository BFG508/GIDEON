function fig_handles = plotQUESTResults(meas, N_STR, STR, DCM_true, DCM_estimated, q_true, q_estimated, quest_info, saveFlag)
%==========================================================================
% plotQUESTResults: Generate comprehensive visualization of QUEST attitude
%                   determination results as separate figures.
%
% Inputs:
%   meas          - Structure array (1xN_STR) with measurements.
%   N_STR         - Number of active star trackers.
%   STR           - STR configuration structure array.
%   DCM_true      - True ECI-to-Body DCM (3x3).
%   DCM_estimated - Estimated ECI-to-Body DCM (3x3).
%   q_true        - True attitude quaternion [qw; qx; qy; qz].
%   q_estimated   - Estimated quaternion [qw; qx; qy; qz].
%   quest_info    - QUEST diagnostic structure from solveQUESTAttitude.
%   saveFlag      - (Optional) Boolean flag to save figures (default: false).
%                   If true, saves to Figures/ folder in FIG, PNG, SVG formats.
%
% Outputs:
%   fig_handles   - Array of figure handles (5 figures).
%
% Generated Figures:
%   1. Star field visualization (STR 1) - pixel coordinates with magnitude
%   2. Residual distribution histogram  - statistical analysis
%   3. DCM error heatmap                - element-wise comparison
%   4. Quaternion comparison            - component-wise bar chart
%   5. Residuals vs. magnitude scatter  - brightness correlation
%==========================================================================

    % Handle optional save_figure input (default: false)
    if nargin < 9
        saveFlag = false;
    end
    
    % Pre-compute common metrics
    q_error = quatmultiply_custom(q_estimated, quatinv_custom(q_true));
    angle_error_arcsec = 2 * acos(min(abs(q_error(1)), 1)) * 206265;
    residuals_arcsec = rad2deg(quest_info.residuals) * 3600;
    
    % Initialize figure handles array
    fig_handles = gobjects(5, 1);
    
    % Create Figures directory if saving is requested
    if saveFlag
        figures_dir = [pwd, '\Figures'];
        if ~exist(figures_dir, 'dir')
            mkdir(figures_dir);
            fprintf('\n✓ Created directory: %s/\n', figures_dir);
        end
        fprintf('\n--- Saving figures ---\n');
    end
    
    %% ====================================================================
    %  FIGURE 1: Star Field Visualization (STR 1)
    % =====================================================================
    
    fig_handles(1) = figure('Name', 'QUEST - Star Field', ...
                            'Color', 'w', 'NumberTitle', 'off');
    
    if N_STR >= 1 && meas(1).N_stars > 0
        % Plot stars
        scatter(meas(1).pixel_x_noisy, meas(1).pixel_y_noisy, ...
                120, meas(1).magnitude, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        
        % Set colormap and limits before colorbar
        colormap(gca, flipud(hot));
        clim([min(meas(1).magnitude), max(meas(1).magnitude)]);
        
        hold on;
        % Add sensor boundaries
        rectangle('Position', [-STR(1).resolution(1)/2, -STR(1).resolution(2)/2, ...
                               STR(1).resolution(1), STR(1).resolution(2)], ...
                  'EdgeColor', 'r', 'LineWidth', 2.5, 'LineStyle', '--');
        
        % Create dummy line for legend (same style as rectangle)
        h_boundary = plot(NaN, NaN, 'r--', 'LineWidth', 2.5);
        hold off;
        
        xlabel('X Pixel', 'FontSize', 13, 'FontWeight', 'bold');
        ylabel('Y Pixel', 'FontSize', 13, 'FontWeight', 'bold');
        title(sprintf('STR 1: %d Stars detected | FOV = %.0f°', ...
                      meas(1).N_stars, STR(1).FOV_deg), ...
              'FontSize', 14, 'FontWeight', 'bold');
        
        grid on;
        axis equal tight;
        xlim([-STR(1).resolution(1)/2, STR(1).resolution(1)/2]);
        ylim([-STR(1).resolution(2)/2, STR(1).resolution(2)/2]);
        
        % Create legend
        legend(h_boundary, 'Sensor Boundary', 'Location', 'northoutside', 'FontSize', 11);
        set(gca, 'FontSize', 11);
        
        % Add colorbar LAST and force rendering
        cb = colorbar;
        cb.Label.String = 'Magnitude [M]';
        cb.Label.FontSize = 12;
        drawnow; % Force MATLAB to complete rendering
    else
        text(0.5, 0.5, 'No Data from STR 1', ...
             'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', [0.5 0.5 0.5]);
        axis off;
    end
    
    if saveFlag
        saveFigure(fig_handles(1), figures_dir, 'QUEST_fig1_star_field');
    end
    
    %% ====================================================================
    %  FIGURE 2: Residual Distribution Histogram
    % =====================================================================
    
    fig_handles(2) = figure('Name', 'QUEST - Residuals', ...
                            'Color', 'w', 'NumberTitle', 'off');
    
    histogram(residuals_arcsec, 30, 'FaceColor', [0.2 0.6 0.9], ...
              'EdgeColor', [0.1 0.3 0.5], 'LineWidth', 1.2, 'FaceAlpha', 0.8);
    
    hold on;
    % Add statistical lines
    xline(mean(residuals_arcsec), 'r-', 'LineWidth', 3, ...
          'Label', sprintf('μ = %.1f"', mean(residuals_arcsec)), ...
          'FontSize', 11, 'LabelVerticalAlignment', 'top', ...
          'LabelHorizontalAlignment', 'left');
    xline(mean(residuals_arcsec) + std(residuals_arcsec), 'r--', 'LineWidth', 2, ...
          'Label', sprintf('μ + σ (%.1f")', mean(residuals_arcsec) + std(residuals_arcsec)), 'FontSize', 10);
    xline(mean(residuals_arcsec) - std(residuals_arcsec), 'r--', 'LineWidth', 2, ...
          'Label', sprintf('μ - σ (%.1f")', mean(residuals_arcsec) - std(residuals_arcsec)), 'FontSize', 10);
    hold off;
    
    xlabel('Residual [arcsec]', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Frequency', 'FontSize', 13, 'FontWeight', 'bold');
    title(sprintf('Residual Distribution | μ = %.1f", σ = %.1f" | N = %d stars', ...
                  mean(residuals_arcsec), std(residuals_arcsec), quest_info.N_obs), ...
          'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 11);
    
    if saveFlag
        saveFigure(fig_handles(2), figures_dir, 'QUEST_fig2_residuals');
    end
    
    %% ====================================================================
    %  FIGURE 3: DCM Error Heatmap
    % =====================================================================
    
    fig_handles(3) = figure('Name', 'QUEST - DCM Error', ...
                            'Color', 'w', 'NumberTitle', 'off');
    
    DCM_error = DCM_estimated - DCM_true;
    DCM_error_scaled = DCM_error * 1e6; % Scale to micro-units
    
    imagesc(DCM_error_scaled);
    
    % Set colormap and limits before colorbar
    colormap(gca, redblueColormap(256));
    clim([-max(abs(DCM_error_scaled(:))), max(abs(DCM_error_scaled(:)))]);
    
    title(sprintf('DCM Error Matrix | Frobenius Norm: %.2e', ...
                  norm(DCM_error, 'fro')), ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Column', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Row', 'FontSize', 13, 'FontWeight', 'bold');
    
    set(gca, 'XTick', 1:3, 'YTick', 1:3, 'FontSize', 12);
    axis square;
    
    % Add text annotations
    for i = 1:3
        for j = 1:3
            if abs(DCM_error_scaled(i,j)) > max(abs(DCM_error_scaled(:))) * 0.5
                text_color = 'w';
            else
                text_color = 'k';
            end
            text(j, i, sprintf('%.2f', DCM_error_scaled(i,j)), ...
                 'HorizontalAlignment', 'center', 'FontSize', 13, ...
                 'Color', text_color, 'FontWeight', 'bold');
        end
    end
    
    % Add colorbar last and force rendering
    cb = colorbar;
    cb.Label.String = 'Error × 10^{-6}';
    cb.Label.FontSize = 12;
    drawnow; % Force MATLAB to complete rendering
    
    if saveFlag
        saveFigure(fig_handles(3), figures_dir, 'QUEST_fig3_DCM_error');
    end
    
    %% ====================================================================
    %  FIGURE 4: Quaternion Comparison
    % =====================================================================
    
    fig_handles(4) = figure('Name', 'QUEST - Quaternion', ...
                            'Color', 'w', 'NumberTitle', 'off');
    
    % Bar plot of quaternion components
    q_error(1) = 1 - q_error(1);
    bar_data = [q_true, q_estimated, q_error];
    b = bar(bar_data, 'grouped', 'BarWidth', 0.85);
    b(1).FaceColor = [0.2 0.7 0.2];
    b(2).FaceColor = [0.2 0.5 0.9];
    b(3).FaceColor = [0.9 0.3 0.2];
    b(1).EdgeColor = 'k';
    b(2).EdgeColor = 'k';
    b(3).EdgeColor = 'k';
    b(1).LineWidth = 1.2;
    b(2).LineWidth = 1.2;
    b(3).LineWidth = 1.2;
    
    xlabel('Quaternion Component', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Value', 'FontSize', 13, 'FontWeight', 'bold');
    title(sprintf('Quaternion Comparison | Angular Error: %.2f arcsec', ...
                  angle_error_arcsec), ...
          'FontSize', 14, 'FontWeight', 'bold');
    legend('True', 'Estimated', 'Error', 'Location', 'northoutside', ...
           'Orientation', 'horizontal', 'FontSize', 12, 'Box', 'on');
    
    set(gca, 'XTickLabel', {'q_w', 'q_x', 'q_y', 'q_z'}, 'FontSize', 12);
    grid on;
    ylim([min(bar_data(:))*1.1, max(bar_data(:))*1.1]);
    
    if saveFlag
        saveFigure(fig_handles(4), figures_dir, 'QUEST_fig4_quaternion');
    end
    
    %% ====================================================================
    %  FIGURE 5: Residuals vs. Star Magnitude
    % =====================================================================
    
    fig_handles(5) = figure('Name', 'QUEST - Residuals vs Magnitude', ...
                            'Color', 'w', 'NumberTitle', 'off');
    
    % Collect all magnitudes
    N_total_stars = sum([meas.N_stars]);
    all_magnitudes = zeros(N_total_stars, 1);
    
    idx_start = 1;
    for i_str = 1:N_STR
        if meas(i_str).N_stars > 0
            idx_end = idx_start + meas(i_str).N_stars - 1;
            all_magnitudes(idx_start:idx_end) = meas(i_str).magnitude;
            idx_start = idx_end + 1;
        end
    end

    % Add trend line with statistics
    p = polyfit(all_magnitudes, residuals_arcsec, 1);
    mag_range = linspace(min(all_magnitudes), max(all_magnitudes), 100);
    h_trend = plot(mag_range, polyval(p, mag_range), 'r--', 'LineWidth', 3);

    hold on;

    % Plot scatter
    scatter(all_magnitudes, residuals_arcsec, 100, residuals_arcsec, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
    
    % Set colormap before colorbar
    colormap(gca, parula);
    
    % Compute correlation
    R = corrcoef(all_magnitudes, residuals_arcsec);
    R_squared = R(1,2)^2;
    
    hold off;
    
    xlabel('Star Magnitude [M]', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Residual [arcsec]', 'FontSize', 13, 'FontWeight', 'bold');
    title('Residuals vs. Star Brightness', 'FontSize', 14, 'FontWeight', 'bold');
    
    grid on;
    set(gca, 'FontSize', 11);
    
    % Create legend with trend line only
    if p(2) >= 0
        legend(h_trend, sprintf('Trend: y = %.2fx + %.2f (R²=%.3f)', p(1), p(2), R_squared), ...
               'Location', 'northoutside', 'FontSize', 11);
    else
        legend(h_trend, sprintf('Trend: y = %.2fx - %.2f (R²=%.3f)', p(1), -p(2), R_squared), ...
               'Location', 'northoutside', 'FontSize', 11);
    end
    
    % Add colorbar last and force rendering
    cb = colorbar;
    cb.Label.String = 'Residual [arcsec]';
    cb.Label.FontSize = 12;
    drawnow; % Force MATLAB to complete rendering
    
    if saveFlag
        saveFigure(fig_handles(5), figures_dir, 'QUEST_fig5_residuals_vs_mag');
    end
    
end