function fig = plotTRIADResults(meas, nSTR, STR, DCMTrue, DCMEst, qTrue, qEst, triadInfo, saveFlag)
%==========================================================================
% plotQUESTResults: Generate comprehensive visualization of QUEST attitude
%                   determination results as separate figures.
%
% Inputs:
%   meas      - Structure array (1xnSTR) with measurements.
%   nSTR      - Number of active star trackers.
%   STR       - STR configuration structure array.
%   DCMTrue   - True ECI-to-Body DCM (3x3).
%   DCMEst    - Estimated ECI-to-Body DCM (3x3).
%   qTrue     - True attitude quaternion [qw; qx; qy; qz].
%   qEst      - Estimated quaternion [qw; qx; qy; qz].
%   triadInfo - TRIAD diagnostic structure from solveTRIADAttitude.
%   saveFlag  - (Optional) Boolean flag to save figures (default: false).
%                If true, saves to Figures/ folder in FIG, PNG, SVG formats.
%
% Outputs:
%   fig       - Array of figure handles (6 figures).
%
% Generated Figures:
%   1. Star field visualization        - pixel coordinates with magnitude
%   2. Residual distribution histogram - statistical analysis
%   3. DCM error heatmap               - element-wise comparison
%   4. Quaternion comparison           - component-wise bar chart
%   5. Residuals vs. magnitude scatter - brightness correlation
%   6. TRIAD vector selection           - primary/secondary vectors and geometry
%==========================================================================

    % Handle optional saveFlag input (default: false)
    if nargin < 9
        saveFlag = false;
    end
    
    % Pre-compute common metrics
    qErr      = quatmultiply(qEst, quatinv(qTrue));
    angErr    = 2 * acos(min(abs(qErr(1)), 1)) * 206265;
    residuals = rad2deg(triadInfo.residuals) * 3600;
    
    % Initialize figure handles array
    fig = gobjects(5, 1);
    
    % Create Figures directory if saving is requested
    if saveFlag
        figDir = fullfile(pwd, 'Figures', 'STR', 'TRIAD');
        if ~exist(figDir, 'dir')
            mkdir(figDir);
            fprintf('\n✓ Created directory: %s/\n', figDir);
        end
        fprintf('\n--- Saving figures ---\n');
    end
    
    %% ====================================================================
    %  FIGURE 1: Star Field Visualization
    % =====================================================================
    
    % Adjust figure size based on number of active STRs
    if nSTR == 1
        fig(1) = figure('Name', 'QUEST - Star Field', ...
                                'Color', 'w', 'NumberTitle', 'off', ...
                                'Position', [100, 100, 700, 600]);
    else
        fig(1) = figure('Name', 'QUEST - Star Field', ...
                                'Color', 'w', 'NumberTitle', 'off', ...
                                'Position', [100, 100, 1400, 600]);
    end
    
    % Determine global magnitude limits for consistent colormap
    magMin =  inf;
    magMax = -inf;
    for iSTR = 1:nSTR
        if meas(iSTR).nStars > 0
            magMin = min([magMin; meas(iSTR).magnitude]);
            magMax = max([magMax; meas(iSTR).magnitude]);
        end
    end
    % Fallback if no stars detected
    if isinf(magMin)
        magMin = 0;
        magMax = 6;
    end
    
    % Dummy handle for shared legend
    hBoundary = gobjects(1);
    
    % --- Loop over active STRs only ---
    for iSTR = 1:nSTR
        % Create subplot only if multiple STRs
        if nSTR > 1
            subplot(1, nSTR, iSTR);
        end
        
        % Check if STR has data
        if meas(iSTR).nStars > 0
            % Plot stars
            scatter(meas(iSTR).pixelX_noisy, meas(iSTR).pixelY_noisy, ...
                    120, meas(iSTR).magnitude, 'filled', ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
            
            % Set colormap and limits
            colormap(gca, flipud(hot));
            clim([magMin, magMax]);
            
            hold on;
            % Add STR boundaries
            rectangle('Position', [-STR(iSTR).resolution(1)/2, -STR(iSTR).resolution(2)/2, ...
                                    STR(iSTR).resolution(1)  ,  STR(iSTR).resolution(2)], ...
                      'EdgeColor', 'r', 'LineWidth', 2.5, 'LineStyle', '--');
            
            % Create dummy line for shared legend (only in first subplot)
            if iSTR == 1
                hBoundary = plot(NaN, NaN, 'r--', 'LineWidth', 2.5);
            end
            hold off;
            
            xlabel('X Pixel', 'FontSize', 13, 'FontWeight', 'bold');
            ylabel('Y Pixel', 'FontSize', 13, 'FontWeight', 'bold');
            
            % Adjust title based on number of STRs
            if nSTR == 1
                title(sprintf('Star Tracker: %d Stars | FOV = %.0f°', ...
                              meas(iSTR).nStars, STR(iSTR).FOV), ...
                      'FontSize', 14, 'FontWeight', 'bold');
            else
                title(sprintf('STR %d: %d Stars | FOV = %.0f°', ...
                              iSTR, meas(iSTR).nStars, STR(iSTR).FOV), ...
                      'FontSize', 14, 'FontWeight', 'bold');
            end
            
            grid on;
            axis equal tight;
            xlim([-STR(iSTR).resolution(1)/2, STR(iSTR).resolution(1)/2]);
            ylim([-STR(iSTR).resolution(2)/2, STR(iSTR).resolution(2)/2]);
            set(gca, 'FontSize', 11);
            
        else
            % Display no data message
            text(0.5, 0.5, sprintf('No Data from STR %d', iSTR), ...
                 'HorizontalAlignment', 'center', 'FontSize', 16, ...
                 'Color', [0.5 0.5 0.5]);
            axis off;
        end
    end
    
    % --- Shared Colorbar ---
    if nSTR == 1
        cb = colorbar('Location', 'eastoutside');
        cb.Label.String = 'Magnitude [M]';
        cb.Label.FontSize = 12;
        cb.Label.FontWeight = 'bold';
    else
        cb = colorbar('Position', [0.92, 0.15, 0.02, 0.7]);
        cb.Label.String = 'Magnitude [M]';
        cb.Label.FontSize = 12;
        cb.Label.FontWeight = 'bold';
    end
    
    % --- Shared Legend ---
    lgd = legend(hBoundary, 'Sensor Boundary', ...
                 'Location', 'northoutside', 'FontSize', 11);
    if nSTR == 2
        lgd.Position = [0.42, 0.92, 0.15, 0.05]; % Centered at top
    end
    drawnow; % Force rendering
    
    if saveFlag
        saveFigure(fig(1), figDir, 'QUEST_fig1_star_field');
    end
    
    %% ====================================================================
    %  FIGURE 2: Residual Distribution Histogram
    % =====================================================================
    
    fig(2) = figure('Name', 'QUEST - Residuals', ...
                            'Color', 'w', 'NumberTitle', 'off');
    
    histogram(residuals, 30, 'FaceColor', [0.2 0.6 0.9], ...
              'EdgeColor', [0.1 0.3 0.5], 'LineWidth', 1.2, 'FaceAlpha', 0.8);
    
    hold on;
    % Add statistical lines
    xline(mean(residuals), 'r-', 'LineWidth', 3, ...
          'Label', sprintf('μ = %.1f"', mean(residuals)), ...
          'FontSize', 11, 'LabelVerticalAlignment', 'top', ...
          'LabelHorizontalAlignment', 'left');
    xline(mean(residuals) + std(residuals), 'r--', 'LineWidth', 2, ...
          'Label', sprintf('μ + σ (%.1f")', mean(residuals) + std(residuals)), 'FontSize', 10);
    xline(mean(residuals) - std(residuals), 'r--', 'LineWidth', 2, ...
          'Label', sprintf('μ - σ (%.1f")', mean(residuals) - std(residuals)), 'FontSize', 10);
    hold off;
    
    xlabel('Residual [arcsec]', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Frequency', 'FontSize', 13, 'FontWeight', 'bold');
    title(sprintf('Residual Distribution | μ = %.1f", σ = %.1f" | N = %d stars', ...
                  mean(residuals), std(residuals), triadInfo.nObs), ...
          'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 11);
    
    if saveFlag
        saveFigure(fig(2), figDir, 'QUEST_fig2_residuals');
    end
    
    %% ====================================================================
    %  FIGURE 3: DCM Error Heatmap
    % =====================================================================
    
    fig(3) = figure('Name', 'QUEST - DCM Error', ...
                            'Color', 'w', 'NumberTitle', 'off');
    
    DCMErr        = DCMEst - DCMTrue;
    DCM_ErrScaled = DCMErr * 1e6; % Scale to micro-units
    
    imagesc(DCM_ErrScaled);
    
    % Set colormap and limits before colorbar
    colormap(gca, redblueColormap(256));
    clim([-max(abs(DCM_ErrScaled(:))), max(abs(DCM_ErrScaled(:)))]);
    
    title(sprintf('DCM Error Matrix | Frobenius Norm: %.2e', ...
                  norm(DCMErr, 'fro')), ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Column', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Row', 'FontSize', 13, 'FontWeight', 'bold');
    
    set(gca, 'XTick', 1:3, 'YTick', 1:3, 'FontSize', 12);
    axis square;
    
    % Add text annotations
    for i = 1:3
        for j = 1:3
            if abs(DCM_ErrScaled(i,j)) > max(abs(DCM_ErrScaled(:))) * 0.5
                text_color = 'w';
            else
                text_color = 'k';
            end
            text(j, i, sprintf('%.2f', DCM_ErrScaled(i,j)), ...
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
        saveFigure(fig(3), figDir, 'QUEST_fig3_DCM_error');
    end
    
    %% ====================================================================
    %  FIGURE 4: Quaternion Comparison
    % =====================================================================
    
    fig(4) = figure('Name', 'QUEST - Quaternion', ...
                            'Color', 'w', 'NumberTitle', 'off');
    
    % Force same sign convention for visualization (use qTrue as reference)
    qEst_plot = qEst;
    if sign(qTrue(1)) ~= sign(qEst(1))
        qEst_plot = -qEst; % Flip sign to match qTrue
    end
    
    % Recompute error with consistent signs
    qErr_plot = quatmultiply(qEst_plot, quatinv(qTrue));

    % Bar plot of quaternion components
    qErr_plot(1)   = 1 - qErr_plot(1);
    barData        = [qTrue, qEst_plot, qErr_plot];
    b              = bar(barData, 'grouped', 'BarWidth', 0.85);
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
                  angErr), ...
          'FontSize', 14, 'FontWeight', 'bold');
    legend('True', 'Estimated', 'Error', 'Location', 'northoutside', ...
           'Orientation', 'horizontal', 'FontSize', 12, 'Box', 'on');
    
    set(gca, 'XTickLabel', {'q_w', 'q_x', 'q_y', 'q_z'}, 'FontSize', 12);
    grid on;
    ylim([min(barData(:))*1.1, max(barData(:))*1.1]);
    
    if saveFlag
        saveFigure(fig(4), figDir, 'QUEST_fig4_quaternion');
    end
    
    %% ====================================================================
    %  FIGURE 5: Residuals vs. Star Magnitude
    % =====================================================================
    
    fig(5) = figure('Name', 'QUEST - Residuals vs Magnitude', ...
                            'Color', 'w', 'NumberTitle', 'off');
    
    % Collect all magnitudes
    nTotalStars = sum([meas.nStars]);
    magnitudes = zeros(nTotalStars, 1);
    
    idxStart = 1;
    for iSTR = 1:nSTR
        if meas(iSTR).nStars > 0
            idxEnd = idxStart + meas(iSTR).nStars - 1;
            magnitudes(idxStart:idxEnd) = meas(iSTR).magnitude;
            idxStart = idxEnd + 1;
        end
    end

    % Add trend line with statistics
    p = polyfit(magnitudes, residuals, 1);
    magRange = linspace(min(magnitudes), max(magnitudes), 100);
    hTrend = plot(magRange, polyval(p, magRange), 'r--', 'LineWidth', 3);
    hold on;

    % Plot scatter
    scatter(magnitudes, residuals, 100, residuals, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
    
    % Set colormap before colorbar
    colormap(gca, parula);
    
    % Compute correlation
    R  = corrcoef(magnitudes, residuals);
    R2 = R(1,2)^2;
    
    hold off;
    
    xlabel('Star Magnitude [M]', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Residual [arcsec]', 'FontSize', 13, 'FontWeight', 'bold');
    title('Residuals vs. Star Brightness', 'FontSize', 14, 'FontWeight', 'bold');
    
    grid on;
    set(gca, 'FontSize', 11);
    
    % Create legend with trend line only
    if p(2) >= 0
        legend(hTrend, sprintf('Trend: y = %.2fx + %.2f (R² = %.3f)', p(1), p(2), R2), ...
               'Location', 'northoutside', 'FontSize', 11);
    else
        legend(hTrend, sprintf('Trend: y = %.2fx - %.2f (R² = %.3f)', p(1), -p(2), R2), ...
               'Location', 'northoutside', 'FontSize', 11);
    end
    
    % Add colorbar last and force rendering
    cb = colorbar;
    cb.Label.String = 'Residual [arcsec]';
    cb.Label.FontSize = 12;
    drawnow; % Force MATLAB to complete rendering
    
    if saveFlag
        saveFigure(fig(5), figDir, 'QUEST_fig5_residuals_vs_mag');
    end
    
    %% ====================================================================
    %  FIGURE 6: TRIAD Vector Selection
    % =====================================================================
    
    fig(6) = figure('Name', 'TRIAD - Vector Selection', ...
                            'Color', 'w', 'NumberTitle', 'off', ...
                            'Position', [200, 100, 1400, 600]);
    
    % Extract indices of primary and secondary vectors
    idxPrimary   = triadInfo.idxPrimary;
    idxSecondary = triadInfo.idxSecondary;
    
    % Reconstruct all observations for visualization
    nTotalObs   = triadInfo.nObs;
    weightsAll  = triadInfo.weights;
    residualsAll = residuals; % Already computed at top of function
    
    % --- SUBPLOT 1: Weight Distribution ---
    subplot(1, 2, 1);
    
    % Plot all observations
    scatter(1:nTotalObs, weightsAll * 100, 120, [0.7 0.7 0.7], 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.0);
    hold on;
    
    % Highlight primary (green) and secondary (blue)
    scatter(idxPrimary, weightsAll(idxPrimary) * 100, 250, [0.2 0.8 0.2], 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 2.5, 'Marker', 'pentagram');
    scatter(idxSecondary, weightsAll(idxSecondary) * 100, 250, [0.2 0.5 0.9], 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 2.5, 'Marker', 'pentagram');
    
    % Add labels
    text(idxPrimary, weightsAll(idxPrimary) * 100 - 1, 'Primary', ...
         'HorizontalAlignment', 'center', 'FontSize', 11, ...
         'FontWeight', 'bold', 'Color', [0.0 0.5 0.0]);
    text(idxSecondary, weightsAll(idxSecondary) * 100 - 1, 'Secondary', ...
         'HorizontalAlignment', 'center', 'FontSize', 11, ...
         'FontWeight', 'bold', 'Color', [0.0 0.3 0.7]);
    
    hold off;
    
    xlabel('Observation Index', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Weight [%]', 'FontSize', 13, 'FontWeight', 'bold');
    title(sprintf('TRIAD Vector Selection | Used: 2/%d observations', nTotalObs), ...
          'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 11);
    
    % --- SUBPLOT 2: Residual Comparison ---
    subplot(1, 2, 2);
    
    barDataResiduals = [residualsAll(idxPrimary); ...
                        residualsAll(idxSecondary); ...
                        mean(residualsAll)];
    bRes = bar(barDataResiduals, 'FaceColor', 'flat');
    bRes.CData(1,:) = [0.2 0.8 0.2];
    bRes.CData(2,:) = [0.2 0.5 0.9];
    bRes.CData(3,:) = [0.7 0.7 0.7];
    bRes.EdgeColor = 'k';
    bRes.LineWidth = 1.2;
    
    ylabel('Residual [arcsec]', 'FontSize', 13, 'FontWeight', 'bold');
    title(sprintf('Residual Analysis | θ = %.1f°', triadInfo.angleBetween), ...
          'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'XTickLabel', {'Primary', 'Secondary', 'All (mean)'}, 'FontSize', 11);
    grid on;
    
    % Add text annotations on bars
    for i = 1:3
        text(i, barDataResiduals(i) + max(barDataResiduals)*0.05, ...
             sprintf('%.2f"', barDataResiduals(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
    end
    
    if saveFlag
        saveFigure(fig(6), figDir, 'TRIAD_fig6_vector_selection');
    end
end