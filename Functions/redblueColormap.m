function cmap = redblueColormap(n)
%==========================================================================
% redblue: Generate diverging red-blue colormap for visualizing positive
%          and negative data (e.g., errors, deviations from zero).
%
% Inputs:
%   n  - (Optional) Number of colormap levels (default: 256).
%        Must be even for symmetric color distribution.
%
% Outputs:
%   cmap - Nx3 colormap array with RGB values in [0, 1].
%          Blue (negative) → White (zero) → Red (positive)
%
% Color Scheme:
%   - Blue end:   RGB = [0, 0, 1] (most negative)
%   - White mid:  RGB = [1, 1, 1] (zero point)
%   - Red end:    RGB = [1, 0, 0] (most positive)
%==========================================================================

    % Set default number of levels if not specified
    if nargin < 1
        n = 256;
    end
    
    % Red channel: ramps from 0 to 1 (left half), then stays at 1 (right half)
    r = [linspace(0, 1, n/2), ones(1, n/2)];
    
    % Green channel: ramps 0→1 (left), then 1→0 (right) for white center
    g = [linspace(0, 1, n/2), linspace(1, 0, n/2)];
    
    % Blue channel: stays at 1 (left half), then ramps 1→0 (right half)
    b = [ones(1, n/2), linspace(1, 0, n/2)];
    
    % Assemble RGB triplets as column vectors (Nx3 array)
    cmap = [r', g', b'];

end
