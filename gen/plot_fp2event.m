function [align,time,events,fig,sp] = plot_fp2event(beh,varargin)
%Align photometry to event times (reward, movement) and plot average trace
%
%   plot_fp2event(beh)
%   [align,time] = plot_fp2event(beh)
%   [align,time,events,fig,sp] = plot_fp2event(beh)
%   [align,time,events,fig,sp] = plot_fp2event(beh,window,flag)
%
%   Description: This function is for aligning photometry signal to
%   discrete events (e.g. reward, movement) for multiple recordings using
%   input structure 'beh'.
%   User will be prompted to respond to two pop-up windows and indicate
%   first which event the photometry signal(s) should be aligned to and
%   second whether average trace should be presented as dF/F% or z-score.
%   The output of this function is figure with nSubplots = nRecordings and 
%   a cell array with the data used to generate each subplot.
%
%   Input:
%   - beh - A data structure containing extracts photometry, behavior
%   information for multiple recordings, using extractBeh or extractData
%   - optional inputs:
%       - window - A vector with window for extraction, in seconds
%       default is [-1 1]
%       - flag - 0 or 1, to indicate whether to plot aligned data
%
%   Output:
%   - align - A cell array with photometry aligned to event times
%   for each recording analyzed
%   - time - A vector of time points for plotting purposes
%   - events - A cell array of event times that photometry is aligned to,
%   in seconds
%   - fig - Figure handle
%   - sp - Subplots handle
%
%   Author: Anya Krok, February 2020

%% Input Variables
evName = menu('Choose Event to Align To:',...
    'Movement Onset','Rest Onset',...
    'Movement Offset','Rest Offset',...
    'Acceleration Peaks','1st Acceleration',...
    'Cue Onset','Reward Delivery','First Lick');
zOrDF = 1; % zOrDF = menu('Choose y-axis:','dF/F %','z-score dF/F');
switch nargin
    case 1
        win = [-1 1]; % Default window is [-1sec : +1sec]
        flag = 1;
    case 2
        if length(varargin{1}) == 2
            win = varargin{1}; flag = 1;
        elseif length(varargin{1}) == 1
            flag = varargin{1};
        end
    case 3
        win = varargin{1};
        flag = varargin{2};
end

%% Align photometry to acceleration times
Fs = beh(1).Fs; % Sampling rate
time = [win(1):1/Fs:win(2)]; % Window for STA
align = cell(length(beh),length(beh(1).FP)); %nShuff = 10; 
events = cell(length(beh),1);
h = waitbar(0, 'STA: signal to acc pks');
for x = 1:length(beh)
    if all(logical(~rem(beh(x).on,1))); diffFs = 1; else; diffFs = 50; end
    switch evName
        case 1; ev = beh(x).on/(Fs/diffFs);     % Extract movement onset times, adjusting event times to be in seconds if necessary
        case 2; ev = beh(x).onRest/(Fs/diffFs); % Extract rest onset times, adjusting event times to be in seconds if necessary
        case 3; ev = beh(x).off/(Fs/diffFs);    % Extract movement offset times, adjusting event times to be in seconds if necessary
        case 4; ev = beh(x).offRest/(Fs/diffFs); % Extract rest offset times, adjusting event times to be in seconds if necessary
        case 5  % Extract acceleration peak times, adjusting event times to be in seconds if necessary
            acc = getAcc(beh(x).vel); % Extract acceleration signal
            [~,locs] = findpeaks(acc,'MinPeakProminence',0.5,'MinPeakDistance',0.5); % Location of peaks, using findpeaks function
            ev = beh(x).time(locs); % Convert peak locations to seconds
            clc         
        case 6 % Extract first acceleration peak times, adjusting event times to be in seconds if necessary
            ev = nan(length(beh(x).on),1);
            for ii = 1:length(beh(x).on) % iterate over movement bouts
                acc = getAcc(beh(x).vel(beh(x).on(ii)*diffFs:beh(x).off(ii)*diffFs)); % extract acceleration during specified movement bout
                [~,locs] = findpeaks(acc,'MinPeakProminence',0.5,'MinPeakDistance',0.5); % location of peaks, using findpeaks function
                if ~isempty(locs)
                    locs_sec = beh(x).time(locs + beh(x).on(ii)*diffFs - 1); % acceleration peaks during movement bout, in seconds
                    ev(ii) = locs_sec(1); % first acceleration in movement bout, in seconds
                end
            end
        case 7; ev = beh(x).cue/(Fs/diffFs);    % Extract cue times, adjusting event times to be in seconds if necessary
        case 8
            if isempty(beh(x).reward); continue; end
            rew = beh(x).reward/(Fs/diffFs); % Extract reward delivery times, adjusting event times to be in seconds
            lick = beh(x).lick/(Fs/diffFs);
            bin = 1/1000; window = [0 1];
            peth = getClusterPETH(lick, rew, bin, window); % PETH: lick aligned to reward in 1 ms bins
            cts = peth.cts{1}; % Lick counts in 1ms bins for each reward trial
            rew_lick0 = find(sum(cts,1) == 0); % Find reward index where total licks within window is 0
            
            bin = 1/1000; window = [-0.25 0];
            peth = getClusterPETH(lick, rew, bin, window);
            rew_prewlick = find(sum(peth.cts{1},1) >= 1); % Find reward index for trials where mouse licks preceding reward
            
            rew([rew_lick0, rew_prewlick]) = nan; 
            cts(:, [rew_lick0, rew_prewlick]) = nan; % Remove non-rewarded trials and trials where mouse licks preceding reward
            ev = rew; % Event time is rewarded trials

        case 9
            rew = beh(x).reward/(Fs/diffFs); % Extract reward delivery times, adjusting event times to be in seconds
            lick = beh(x).lick/(Fs/diffFs);
            bin = 1/1000; window = [0 1];
            peth = getClusterPETH(lick, rew, bin, window); % PETH: lick aligned to reward in 1 ms bins
            cts = peth.cts{1}; % Lick counts in 1ms bins for each reward trial
            rew_lick0 = find(sum(cts,1) == 0); % Find reward index where total licks within window is 0
            rew(rew_lick0) = nan; cts(:, rew_lick0) = nan; % Remove non-rewarded trials
            [~, firstLick] = max(cts~=0, [], 1); % Find first non-zero index for each trial
            ev = rew(:) + firstLick(:).*bin; % Event time is first lick after reward for rewarded trials only
            % ev = lick;   % Extract lick times, adjusting event times to be in seconds if necessary
    end
    events{x} = ev;
    if isempty(beh(x).FP); continue; end %if no photometry, continue to next recording
    for y = 1:length(beh(x).FP)
        sig = beh(x).FP{y}; % Signal that will be aligned to event times
        [mat,~,mat_z] = getSTA(sig, ev, Fs, [time(1), time(end)]); % STA: aligning photometry to event times
%         ev_new = shiftST(ev, nShuff, 1/nShuff); % Shuffle event times n times
%         mat_shuff = [];
%         for z = 1:nShuff
%             tmp_shuff = getSTA(sig, ev_new{z}, Fs, [time(1), time(end)]); % STA: aligning photometry to shuffled event times
%             mat_shuff = [mat_shuff, tmp_shuff]; 
%         end
%         mu = nanmean(nanmean(mat_shuff,2)); sigma = nanmean(nanstd(mat_shuff,[],2)); % Use mu, sigma of shuffled STA to z-score 
        switch zOrDF
            case 1
                align{x,y} = mat; ylbl = 'dF/F %'; % Load STA into output cell array
            case 2
                align{x,y} = mat_z; ylbl = 'z-score'; % Load STA into output cell array
        end
    end
    waitbar(x/length(beh),h);
end
close(h); fprintf('Done aligning photometry to events! \n');

%% Plot for each recording
if flag == 1
    fig = figure; % Figure handle
    plm = floor(sqrt(size(align,1))); pln = ceil(size(align,1)/plm); % Subplot size depending on number of recordings
    clr = {'g','m','b','r'};
    switch evName % Assign label based on event time aligning to
        case 1; lbl = 'Movement Onset'; case 2; lbl = 'Rest Onset'; 
        case 3; lbl = 'Movement Offset'; case 4; lbl = 'Rest Offset';
        case 5; lbl = 'Acceleration Peak'; case 6; lbl = '1st Acceleration';
        case 7; lbl = 'Cue Onset'; case 8; lbl = 'Reward Delivery'; case 9; lbl = 'Lick';
    end
    for x = 1:size(align,1) % Iterate over each recording
        sp(x) = subplot(plm,pln,x); 
        for y = 1:size(align,2) % Iterate over each photometry signal
            if isempty(align{x,y}); continue; end
            shadederrbar(time, nanmean(align{x,y},2), SEM(align{x,y},2), clr{y}); hold on % Plot average STA for each photometry signal
        end
        xlabel(sprintf('Latency to %s (s)',lbl)); 
        ylabel(sprintf('FP (%s)',ylbl)); grid on; xlim([time(1) time(end)]);
        title(sprintf('%s - %s',beh(x).rec,beh(x).site)); 
    end
    % linkaxes(sp,'y'); 
    linkaxes(sp,'x'); % Link x,y axes
end
end
