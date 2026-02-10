function [taus, sigma_allan] = computeAllanDeviation(x, dt)
%==========================================================================
% computeAllanDeviation: Calculate the Allan Deviation (ADEV) for a given
%                        time series data.
%
% Description:
%   The Allan Deviation is a two-sample variance metric used to characterize
%   noise stability in clocks, oscillators, and inertial sensors (gyros/accels).
%   It helps distinguish between different noise types (white noise, random 
%   walk, flicker noise) based on the slope of the log-log plot.
%
% Inputs:
%   x           - Input time series data (1xN or Nx1 vector).
%                 (Ideally zero-mean error signals like gyro rate error).
%   dt          - Sampling interval [s].
%
% Outputs:
%   taus        - Averaging time intervals [s].
%   sigma_allan - Allan Deviation values corresponding to each tau.
%
% Algorithm:
%   1. Define a set of averaging cluster sizes (m).
%   2. For each cluster size m:
%      a. Divide data into K adjacent, non-overlapping blocks of length m.
%      b. Compute the mean value of each block.
%      c. Calculate the differences between consecutive block means.
%      d. Compute the variance of these differences (Allan Variance).
%   3. Take the square root to get Allan Deviation.
%==========================================================================

    % Ensure input is a column vector
    x = x(:);
    N = length(x);
    
    % Determine maximum cluster size (m)
    % ADEV becomes statistically unreliable for small K (number of clusters).
    % We limit the max averaging time to N/10 to ensure at least ~10 samples.
    max_m = floor(N/10);
    
    % Edge case: insufficient data points for analysis
    if max_m < 2
        warning('Insufficient data points for Allan Deviation. Returning standard deviation.');
        taus = dt;
        sigma_allan = std(x);
        return;
    end
    
    % Generate logarithmically spaced cluster sizes (m)
    % This provides uniform density on a log-log plot.
    % We generate ~30 points spanning from 1 sample to max_m.
    m_list = unique(round(logspace(0, log10(max_m), 30)));
    
    % Initialize outputs
    taus = m_list * dt;       % Convert samples to time [s]
    sigma_allan = zeros(size(taus));
    
    % Loop over each averaging cluster size
    for k = 1:length(m_list)
        m = m_list(k);
        
        % Number of non-overlapping clusters of length m
        K = floor(N / m);
        
        % Safety check (should be covered by max_m logic, but robust)
        if K < 2
            sigma_allan(k) = NaN;
            continue;
        end
        
        % Compute cluster averages (vectorized for speed could be done via reshape)
        % Note: Using loop here for clarity, but 'reshape' + 'mean' is faster for perfect multiples.
        % The manual loop handles the truncation of the trailing end of data naturally.
        y = zeros(K,1);
        for i = 1:K
            idx_start = (i-1)*m + 1;
            idx_end   = i*m;
            y(i) = mean(x(idx_start:idx_end));
        end
        
        % Calculate Allan Variance: sigma^2 = 0.5 * E[ (y_{i+1} - y_i)^2 ]
        dy = diff(y);                 % Differences between consecutive averages
        sigma2 = 0.5 * mean(dy.^2);   % Mean squared difference / 2
        
        % Allan Deviation is the square root
        sigma_allan(k) = sqrt(sigma2);
    end
    
    % Filter out any NaNs resulting from edge cases
    valid = ~isnan(sigma_allan);
    taus = taus(valid);
    sigma_allan = sigma_allan(valid);

end