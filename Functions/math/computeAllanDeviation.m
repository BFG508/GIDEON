function [tau, sigmaAllan] = computeAllanDeviation(x, dt)
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
%   x          - Input time series data (1xN or Nx1 vector).
%                (Ideally zero-mean error signals like gyro rate error).
%   dt         - Sampling interval [s].
%
% Outputs:
%   tau        - Averaging time intervals [s].
%   sigmaAllan - Allan Deviation values corresponding to each tau.
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
    % The max averaging time is limit to N/10 to ensure at least ~ 10 samples.
    maxM = floor(N/10);
    
    % Edge case: insufficient data points for analysis
    if maxM < 2
        warning('Insufficient data points for Allan Deviation. Returning standard deviation.');
        tau        = dt;
        sigmaAllan = std(x);
        return;
    end
    
    % Generate logarithmically spaced cluster sizes (m)
    % This provides uniform density on a log-log plot.
    % ~30 points are generating spanning from 1 sample to maxM.
    mList = unique(round(logspace(0, log10(maxM), 30)));
    
    % Initialize outputs
    tau = mList * dt; % Convert samples to time [s]
    sigmaAllan = zeros(size(tau));
    
    % Loop over each averaging cluster size
    for k = 1:length(mList)
        m = mList(k);
        
        % Number of non-overlapping clusters of length m
        K = floor(N / m);
        
        % Safety check (should be covered by maxM logic, but robust)
        if K < 2
            sigmaAllan(k) = NaN;
            continue;
        end
        
        % Compute cluster averages
        y = mean(reshape(x(1:K*m), m, K), 1);

        % Calculate Allan Variance: sigma^2 = 1/2 * E[ (y_{i+1} - y_i)^2 ]
        dy = diff(y);               % Differences between consecutive averages
        sigma2 = 1/2 * mean(dy.^2); % Mean squared difference / 2
        
        % Allan Deviation is the square root
        sigmaAllan(k) = sqrt(sigma2);
    end
    
    % Filter out any NaNs resulting from edge cases
    valid      = ~isnan(sigmaAllan);
    tau        = tau(valid);
    sigmaAllan = sigmaAllan(valid);
end