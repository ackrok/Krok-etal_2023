function [align,time,events,fig,sp] = plot_fp2fp(beh,varargin)
%Align one photometry signal to peaks of another photometry signal and plot average trace
%
%   plot_fp2fp(beh)
%   [align,time,events] = plot_fp2fp(beh)
%   [align,time,events,fig,sp] = plot_fp2fp(beh)
%   [align,time,events,fig,sp] = plot_fp2fp(beh,window)
%
%   Description: This function is for aligning one photometry signal to
%   peaks of another photometry signal using input structure 'beh'.
%   User will be prompted to respond to two pop-up windows and indicate
%   first which event the photometry peaks(s) should be stratified for and
%   second whether average trace should be presented as dF/F% or z-score.
%   The output of this function is figure with nSubplots = nRecordings and 
%   a cell array with the data used to generate each subplot.
%
%   For ACh, DA data where FP{1} is ACh signal and FP{2} is DA signal, then
%   align{x,1} is DA signal aligned to ACh peaks
%   align{x,2} is ACh signal aligned to DA peaks
%
%   Input:
%   - beh - A data structure containing extracts photometry, behavior
%   information for multiple recordings, using extractBeh or extractData
%   - optional inputs:
%       - window - A vector with window for extraction, in seconds
%       default is [-1 1]
%
%   Output:
%   - align - A cell array with photometry aligned to fp peak times
%   for each recording analyzed
%   - time - A vector of time points for plotting purposes
%   - events - A cell array of event times that photometry is aligned to,
%   in seconds
%   - fig - Figure handle
%   - sp - Subplots handle
%
%   Author: Anya Krok, June 2021

%% Input Variables
evStratify = menu('Choose Event to Stratify By:','All Peaks','Rest','Locomotion','Reward','Lick');
zOrDF = menu('Choose y-axis:','dF/F %','z-score dF/F');
if nargin == 2
    win = varargin{1};
else
    win = [-1 1]; % Default window is [-1sec : +1sec]
end
winRew = [0.2 0.8]; % Window around reward delivery
winLick = [0 0.01]; % Window around lick onset

%% Align photometry to acceleration times
Fs = beh(1).Fs; % Sampling rate
time = [win(1):1/Fs:win(2)]; % Window for STA
align = cell(length(beh),length(beh(1).FP)); %nShuff = 10; 
events = align; % Initialize empty cell array
h = waitbar(0, 'STA: FP signal to FP peaks');
for x = 1:length(beh)
    % if all(logical(~rem(beh(x).on,1))); diffFs = 1; else; diffFs = 50; end
    if isempty(beh(x).FP); continue; end %if no photometry, continue to next recording
    if length(beh(x).FP) < 2; continue; end
    for y = 1:length(beh(x).FP)
        fp = beh(x).FP{y}; % Signal that will be aligned to event times
        peakProm = 0.1; peakDist = 0.5; % Parameters for FP peaks
        fp_norm = (fp - min(fp)) / (max(fp) - min(fp)); % Min Max Normalization
        if strcmp(strtok(beh(x).FPnames{y},'-'),'ACh'); fp_norm = -1*fp_norm; end % Use troughs for ACh
        [~,locs_samp] = findpeaks(fp_norm,'MinPeakProminence',peakProm,'MinPeakDistance',peakDist); 
        locs = beh(x).time(locs_samp);
        pks = fp(locs_samp);
        switch evStratify
            case 1; ev = locs; % All FP peaks
            case 2; ev = extractEventST(locs, beh(x).onRest/Fs, beh(x).offRest/Fs, 1); % FP peak times during rest
            case 3; ev = extractEventST(locs, beh(x).on/Fs, beh(x).off/Fs, 1); % FP peak times during movement
            case 4; ev = extractEventST(locs, beh(x).reward/Fs + winRew(1), beh(x).reward/Fs + winRew(2), 1); % FP peak times during reward
            case 5; ev = extractEventST(locs, beh(x).lick/Fs + winLick(1), beh(x).lick/Fs + winLick(2), 1); % FP peak times during licking
        end
        if isempty(ev); continue; end % If no peaks, continue to next iteration
        switch y
            case 1; sig = beh(x).FP{2}; 
            case 2; sig = beh(x).FP{1}; 
        end
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
        events{x,y} = ev;
    end
    waitbar(x/length(beh),h);
end
close(h); fprintf('Done aligning photometry to events! \n');

%% Plot for each recording
fig = figure; % Figure handle
plm = floor(sqrt(size(align,1))); pln = ceil(size(align,1)/plm); % Subplot size depending on number of recordings
clr = {'m','g','b','r'};
switch evStratify % Assign label based on event time aligning to
    case 1; lbl = 'All Peaks'; 
    case 2; lbl = 'Rest'; case 3; lbl = 'Movement'; 
    case 4; lbl = 'Reward'; case 5; lbl = 'Lick';
end
for x = 1:size(align,1) % Iterate over each recording
    if length(beh(x).FP) < 2; continue; end
    sp(x) = subplot(plm,pln,x); 
    for y = 1:size(align,2) % Iterate over each photometry signal
        if isempty(align{x,y}); continue; end
        shadederrbar(time, nanmean(align{x,y},2), SEM(align{x,y},2), clr{y}); hold on % Plot average STA for each photometry signal
    end
    xlabel(sprintf('Latency to %s peaks (s)',beh(x).FPnames{y})); 
    ylabel(sprintf('FP (%s)',ylbl)); grid on; xlim([time(1) time(end)]);
    title(sprintf('%s - %s - %s',beh(x).rec,beh(x).site,lbl)); 
end
% linkaxes(sp,'y'); 
linkaxes(sp,'x'); % Link x,y axes

end
