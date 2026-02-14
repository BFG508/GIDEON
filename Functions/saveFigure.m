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
    figFile = fullfile(figDir, [figName, '.fig']);
    savefig(figHandle, figFile);
    fprintf('  ✓ Saved: %s\n', figFile);
    
    % Save as PNG raster image (300 DPI for publication quality)
    pngFile = fullfile(figDir, [figName, '.png']);
    print(figHandle, pngFile, '-dpng', '-r300');
    fprintf('  ✓ Saved: %s\n', pngFile);
    
    % Save as SVG vector graphic (scalable, ideal for presentations/LaTeX)
    svgFile = fullfile(figDir, [figName, '.svg']);
    print(figHandle, svgFile, '-dsvg', '-vector');
    fprintf('  ✓ Saved: %s\n', svgFile);

end
