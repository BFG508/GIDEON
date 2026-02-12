function fig = plotMAGResults(t, BTrue_body, meas, MAG, saveFlag)
%==========================================================================
% plotMAGResults: Generate quicklook figures for MAG truth vs measurements.
% Inputs:
%   t             - Time vector                       [s], 1xN or Nx1
%   BTrue_body    - True magnetic field in body frame [nT], 3xN
%   meas          - MAG measurement struct:
%                     .B                              [nT], 3xN
%                     .BClean                         [nT], 3xN
%                     .BDet                           [nT], 3xN
%                     .biasDyn                        [nT], 3xN
%                     .biasTotal                      [nT], 3xN
%   MAG           - MAG parameter struct (optional fields used):
%                     .DCM_mounting                   (3x3)
%                     .resolution                     [nT]
%                     .range                          [nT]
%   saveFlag      - (Optional) If true, saves figures to Figures/MAG/
%                   Default: false
%
% Outputs:
%   fig           - Array of figure handles
%
% Plots generated:
%   1) Magnetic field error (components + norm) in body frame
%   2) Magnetic field direction errorand magnitude error
%   3) Dynamic bias evolution (Gauss-Markov drift)
%   4) 3D B-vector cloud (measured vs truth)
%==========================================================================

    if nargin < 5
        saveFlag = false;
    end
    
    % Ensure time is row vector
    t = t(:).';
    N = numel(t);
    
    fig = gobjects(4,1);
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

    % Frame harmonization
    % Truth in body frame is given: BTrue_body (3xN)
    if ~isequal(size(BTrue_body), [3, N])
        error('BTrue_body must be 3xN with N = length(t).');
    end

    % Measured B is in body frame already (meas.B is final digital output)
    BMeas = meas.B;
    BTrue = BTrue_body;
    
    %% ====================================================================
    % FIGURE 1: Magnetic Field Vector Error (Components + Norm)
    % =====================================================================
    
    Nfig = Nfig + 1;
    fig(Nfig) = figure('Name', 'MAG - B Vector Error', ...
                        'Color', 'w', 'NumberTitle', 'off', ...
                        'Position', [500, 100, 1000, 800]);

    BErr     = BMeas - BTrue;
    BErrNorm = vecnorm(BErr, 2, 1);
    
    axesLabels = {'x','y','z'};
    for i = 1:3
        subplot(4,1,i);
        plot(t, BErr(i,:), 'b-', 'LineWidth', 1.2);
        ylabel(sprintf('\\DeltaB_%s [nT]', axesLabels{i}), 'FontSize', 11, 'FontWeight', 'bold');
        if i == 1
            title('Magnetic Field Error (Body Frame)', 'FontSize', 13, 'FontWeight', 'bold');
        end

        grid on;
        set(gca, 'FontSize', 10);
        xlim([t(1), t(end)]);
    end
    
    subplot(4,1,4);
    plot(t, BErrNorm, 'k-', 'LineWidth', 1.2);
    xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('||\DeltaB|| [nT]', 'FontSize', 11, 'FontWeight', 'bold');

    grid on;
    set(gca, 'FontSize', 10);
    xlim([t(1), t(end)]);
    
    if saveFlag
        saveFigure(fig(Nfig), saveDir, 'MAG_fig1_B_error');
    end
    
    %% ====================================================================
    % FIGURE 2: Direction and Magnitude Error
    % =====================================================================
    
    Nfig = Nfig + 1;
    fig(Nfig) = figure('Name', 'MAG - Direction & Magnitude Error', ...
                       'Color', 'w', 'NumberTitle', 'off');

    uTrue = BTrue ./ vecnorm(BTrue, 2, 1);
    uMeas = BMeas ./ vecnorm(BMeas, 2, 1);
    
    cosAng = max(min(sum(uTrue .* uMeas, 1), 1), -1);
    angErr = rad2deg(acos(cosAng));
    
    BMagTrue = vecnorm(BTrue, 2, 1);
    BMagMeas = vecnorm(BMeas, 2, 1);
    magErr   = BMagMeas - BMagTrue;
    
    subplot(2,1,1);
    plot(t, angErr, 'm-', 'LineWidth', 1.2);
    ylabel('Angle [deg]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Field Direction Error', 'FontSize', 13, 'FontWeight', 'bold');

    set(gca, 'FontSize', 10);
    grid on;
    xlim([t(1), t(end)]);
    
    subplot(2,1,2);
    plot(t, magErr, 'b-', 'LineWidth', 1.2);
    xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('\\Delta||B|| [nT]', 'FontSize', 11, 'FontWeight', 'bold');

    set(gca, 'FontSize', 10);
    grid on;
    xlim([t(1), t(end)]);
    
    if saveFlag
        saveFigure(fig(Nfig), saveDir, 'MAG_fig2_direction_magnitude_error');
    end
    
    %% ====================================================================
    % FIGURE 3: Dynamic Bias Evolution
    % =====================================================================

    Nfig = Nfig + 1;
    fig(Nfig) = figure('Name', 'MAG - Bias Drift', ...
                       'Color', 'w', 'NumberTitle', 'off');
    
    % Bias is in MAG frame, rotate to body for plotting
    biasDynBody = MAG.DCM_mounting' * meas.biasDyn;
    
    hold on;
    plot(t, biasDynBody(1,:), 'r', 'LineWidth', 1.2, 'DisplayName', 'x');
    plot(t, biasDynBody(2,:), 'g', 'LineWidth', 1.2, 'DisplayName', 'y');
    plot(t, biasDynBody(3,:), 'b', 'LineWidth', 1.2, 'DisplayName', 'z');
    hold off;
    
    xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Dynamic Bias [nT]', 'FontSize', 11, 'FontWeight', 'bold');
    title('MAG Dynamic Bias (Gauss-Markov Process)', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northoutside', 'Orientation', 'horizontal', ...
           'FontSize', 10);

    set(gca, 'FontSize', 10);
    grid on;
    xlim([t(1), t(end)]);
    
    if saveFlag
        saveFigure(fig(Nfig), saveDir, 'MAG_fig3_bias_drift');
    end
    
    %% ====================================================================
    % FIGURE 4: 3D B-Vector Cloud + 2D Projections (Measured vs Truth)
    % =====================================================================

    Nfig = Nfig + 1;
    fig(Nfig) = figure('Name', 'MAG - B Cloud (3D + Projections)', ...
                       'Color', 'w', 'NumberTitle', 'off', ...
                       'Position', [100, 100, 1000, 800]);
    
    % Downsample for plotting performance if needed
    idx = 1:N;
    if N > 8000
        idx = round(linspace(1, N, 8000));
    end
    
    % ---- (1) Top-Left: x vs y ----
    subplot(2,2,1);
    plot(BTrue(1,idx), BTrue(2,idx), 'b-', 'MarkerSize', 4, 'DisplayName', 'Truth');
    hold on;
    plot(BMeas(1,idx), BMeas(2,idx), 'r-', 'MarkerSize', 4, 'DisplayName', 'Measured');
    hold off;
    

    xlabel('B_x [nT]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('B_y [nT]', 'FontSize', 11, 'FontWeight', 'bold');
    title('B_x vs B_y', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 10);
    
    % ---- (2) Top-Right: x vs z ----
    subplot(2,2,2);
    plot(BTrue(1,idx), BTrue(3,idx), 'b-', 'MarkerSize', 4, 'DisplayName', 'Truth'); 
    hold on;
    plot(BMeas(1,idx), BMeas(3,idx), 'r-', 'MarkerSize', 4, 'DisplayName', 'Measured'); 
    hold off;

    xlabel('B_x [nT]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('B_z [nT]', 'FontSize', 11, 'FontWeight', 'bold');
    title('B_x vs B_z', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 10);
    
    % ---- (3) Bottom-Left: y vs z ----
    subplot(2,2,3);
    plot(BTrue(2,idx), BTrue(3,idx), 'b-', 'MarkerSize', 4, 'DisplayName', 'Truth');
    hold on;
    plot(BMeas(2,idx), BMeas(3,idx), 'r-', 'MarkerSize', 4, 'DisplayName', 'Measured'); 
    hold off;
    
    xlabel('B_y [nT]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('B_z [nT]', 'FontSize', 11, 'FontWeight', 'bold');
    title('B_y vs B_z', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 10);
    
    % ---- (4) Bottom-Right: 3D cloud ----
    subplot(2,2,4);
    plot3(BTrue(1,idx), BTrue(2,idx), BTrue(3,idx), ...
          'b-', 'MarkerSize', 4, 'DisplayName', 'Truth');
    hold on;
    plot3(BMeas(1,idx), BMeas(2,idx), BMeas(3,idx), ...
          'r-', 'MarkerSize', 4, 'DisplayName', 'Measured'); hold off;

    xlabel('B_x [nT]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('B_y [nT]', 'FontSize', 11, 'FontWeight', 'bold');
    zlabel('B_z [nT]', 'FontSize', 11, 'FontWeight', 'bold');
    title('3D Magnetic Field Cloud', 'FontSize', 12, 'FontWeight', 'bold');

    grid on;
    set(gca, 'FontSize', 10);
    view(35, 25);
    
    % One legend for the whole figure (simple approach: put it in the last axes)
    lgd = legend('Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 10);
    lgd.Position(1) = 0.5 - lgd.Position(3)/2;
    lgd.Position(2) = 0.96;
    
    if saveFlag
        saveFigure(fig(Nfig), saveDir, 'MAG_fig4_3D_cloud');
    end
    
end
