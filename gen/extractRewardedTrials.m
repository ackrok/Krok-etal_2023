function [rewYes, rewNo, lickNew] = extractRewardedTrials(rewDelivery, lick, varargin)
%Extract index for rewarded and non-rewarded trials based on whether animal
%licked within specified window
%
%   [rewYes, rewNo, lickNew] = extractRewardedTrials(rewDelivery, lick)
%   [rewYes, rewNo, lickNew] = extractRewardedTrials(rewDelivery, lick, window)
%
%   Description: This function is for extraction of indices corresponding
%   to trials where reward was and was not collected, according to whether
%   animal licked within specified window
%
%   Input:
%   - rewDelivery - Vector of reward delivery times, in seconds
%   - lick - Vector of lick times, in seconds
%   - window (optional) - Window for reward collection, in seconds
%       default is [0 1] second
%
%   Output:
%   - rewYes - Vector with indices for delivery times where reward was
%   collected
%   - rewNo - Vector with indices for non-rewarded trials
%   - lickNew - Vector of lick times
%
%   Author: Anya Krok, August 2022
    %%
    switch nargin
        case 2
            window = [0 1]; % default window is 1 second
        case 3
            window = varargin{1};
    end
    lick = lick(:);
    rewDelivery = rewDelivery(:);
    
    %% Correct lick event times 
    %to ensure behaviorally possible separation between events
    lick_repeat = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lickNew = [lick(1); lick_sub(lick_repeat)]; % Overwrite lick vector

    %% 
    bin = 1/1000;
    peth = getClusterPETH(lickNew, rewDelivery, bin, window); % PETH: lick aligned to reward in 1 ms bins
    cts = peth.cts{1}; % Lick counts in 1ms bins for each reward trial
    rew_lick0 = find(sum(cts,1) == 0); % Find reward index where total licks within window is 0

    rew_prelick = [];
    bin = 1/1000; window = [-0.2 0];
    peth = getClusterPETH(lickNew, rewDelivery, bin, window);
    rew_prelick = find(sum(peth.cts{1},1) >= 1); % Find reward index for trials where mouse licks preceding reward
        
    rewDelivery([rew_lick0, rew_prelick]) = nan; 
    cts(:, [rew_lick0, rew_prelick]) = nan; % Remove non-rewarded trials and trials where mouse licks preceding reward

    rewNo = find(isnan(rewDelivery)); % Index of deliveries where animal did not lick to receive reward
    rewYes = find(~isnan(rewDelivery)); % Index of delivieries where animal collected reward
end
