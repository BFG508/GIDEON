function saveFigure(figHandle, figDir, figName)
%==========================================================================
% saveFigure: Save MATLAB figure in multiple formats for publication 
%             and archival purposes.
%
% Inputs:
%   figHandle - Handle to MATLAB figure object to be saved.
%   figDir    - Target directory path (string or char array).
%   figName   - Base filename without extension (string or char array).
%
% Outputs:
%   None (generates 3 files on disk: .fig, .png, .svg)
%
% File Formats:
%   - .fig: Native MATLAB format (editable, preserves all properties)
%   - .png: Raster format at 300 DPI (high-quality publications)
%   - .svg: Vector format (scalable, ideal for LaTeX/presentations)
%==========================================================================

    % Save as MATLAB .fig (native format, fully editable)
    fig_file = fullfile(figDir, [figName, '.fig']);
    savefig(figHandle, fig_file);
    fprintf('  ✓ Saved: %s\n', fig_file);
    
    % Save as PNG raster image (300 DPI for publication quality)
    png_file = fullfile(figDir, [figName, '.png']);
    print(figHandle, png_file, '-dpng', '-r300');
    fprintf('  ✓ Saved: %s\n', png_file);
    
    % Save as SVG vector graphic (scalable, ideal for presentations/LaTeX)
    svg_file = fullfile(figDir, [figName, '.svg']);
    print(figHandle, svg_file, '-dsvg', '-vector');
    fprintf('  ✓ Saved: %s\n', svg_file);

end
