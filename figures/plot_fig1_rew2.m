behLoaded = menu('beh loaded into workspace already?','yes','no');

switch behLoaded
    case 2
        %% Select .mat files you want to add to summary data structu
        fPath = 'R:\tritsn01labspace\'; 
        [fName,fPath] = uigetfile([fPath,'*.mat'],'MultiSelect','On');
        beh = extractBeh(fPath, fName);
end
if ~exist('beh','var'); error('No variable called beh exists'); end

%% Variables
lickWithin = 0.25; %CHANGE, lick within this window
winRew = [-1 2]; % CHANGE, window for aligning signal to rewarded trials
winBase = [-4 -1]; % window for baseline
NumStd = 2; % CHANGE, for immobility trough/peak analysis

%% Align to reward
a_beh = cell(length(beh),2); % initialize output
nFP = length(beh(x).FP); % number of photometry signals
for x = 1:length(beh) % iterate over recordings
    Fs = beh(x).Fs;
    rew = beh(x).reward./Fs; % reward delivery times, in seconds
    lick = beh(x).lick./Fs; % lick times, in seconds
    rewYes = extractRewardedTrials(rew, lick, [0 lickWithin]); 
    ev = rew(rewYes); % aliging signal to rewarded trails
    for y = 1:nFP
        sig = beh(x).FP{y}; % signal that will be aligned to event times
        sig = sig - nanmean(sig); % subtract mean of trace to center
        [mat, time] = getSTA(sig, ev, Fs, [winBase(1), winRew(2)]);
        matBase = nanmean(mat(find(time == winBase(1)):find(time == winBase(2)),:),1); % average across baseline window
        a_beh{x,y} = mat - matBase;
    end
end
% 
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end % extract mouse ID
uni = unique(rec);
nAn = length(uni); % number of unique mouse recordings
% 
a_rew = cell(nAn,2); % initialize output
a_rew_avg = cell(1,2); % initialize output
a_rew_norm = cell(1,2); % initialize output
for x = 1:nAn
    ii = find(strcmp(rec,uni{x})); % match mouse ID
    for y = 1:nFP
        a_rew{x,y} = [a_beh{ii,y}]; % concatenate aligned signals from multiple recordings from same animal
        a_rew{x,y} = a_rew{x,y}(find(time == winRew(1)):find(time == winRew(2)),:); % extract reward window
        a_rew_avg{y}(:,x) = nanmean(a_rew{x,y},2); % average across all rewarded trials
    end
    y = 1; tmp = nanmean([a_beh{ii,y}],2);
    bb = nanmean(tmp(find(time == winBase(1)):find(time == winBase(2)),:),1); % baseline
    mm = min(tmp(find(time == 0):find(time == 1),:)); % minimum
    norm = (tmp - mm)./(bb - mm); % normalize
    norm = norm(find(time == winRew(1)):find(time == winRew(2)),:); % extract reward window
    norm = norm - 1; % shift downward such that baseline is at 0 
    a_rew_norm{y}(:,x) = norm;
    
    y = 2; tmp = nanmean([a_beh{ii,y}],2);
    bb = nanmean(tmp(find(time == winBase(1)):find(time == winBase(2)),:),1); % baseline
    mm = max(tmp(find(time == 0):find(time == 1),:)); % minimum
    norm = (tmp - bb)./(mm - bb); % normalize
    norm = norm(find(time == winRew(1)):find(time == winRew(2)),:); % extract reward window
    a_rew_norm{y}(:,x) = norm;
end
t_rew = time(find(time == winRew(1)):find(time == winRew(2))); % extract reward window
m = max(a_rew_avg{2}); [~,sortIdx] = sort(m); % sort by maximal DA peak
% 
%% Extract peak and trough values
lag = []; val = []; % initialize output
for x = 1:size(a_rew,1)
    y = 1; r = [0.1 0.3]; Fs = 50; r = (r/(1/Fs) + find(t_rew==0)); % range for ACh peak
    [a,b] = max(a_rew{x,y}(r(1):r(2),:)); % find local maximum within range
    c = t_rew(b + r(1) - 1); % convert index to seconds
    c(b == 1) = nan; c(isnan(a)) = nan; % remove values that are not local maxima
    lag(x,1) = nanmean(c); % average across rewarded trials for mouse
    val(x,1) = nanmean(a); % average across rewarded trials for mouse
    % 
    y = 2; r = [0.1 0.5]; r = (r/(1/Fs) + find(t_rew==0)); % range for DA peak
    [a,b] = max(a_rew{x,y}(r(1):r(2),:)); % find local maximum within range
    c = t_rew(b + r(1) - 1); c(isnan(a)) = nan; % convert index to seconds
    lag(x,2) = nanmean(c); % average across rewarded trials for mouse
    val(x,2) = nanmean(a); % average across rewarded trials for mouse
    % 
    y = 1; r = [0.2 0.7]; r = (r/(1/Fs) + find(t_rew==0)); % range for ACh trough
    [a,b] = min(a_rew{x,y}(r(1):r(2),:)); % find local minima within range
    c = t_rew(b + r(1) - 1); c(isnan(a)) = nan; % convert index to seconds
    lag(x,3) = nanmean(c); % average across rewarded trials for mouse
    val(x,3) = nanmean(a); % average across rewarded trials for mouse
end
% 
%% PLOTTING
fig = figure; fig.Position([3 4]) = [1375 800];
clr = {'g','m'};
y_lims = [-4 4; -2 10];
for y = 1:2
    plotme = a_rew_avg{y}; 
    fp_lbl = beh(1).FPnames{y};
    subplot(2, 3, 1 + ((y-1)*3)); hold on
        h = imagesc(t_rew, [1:size(plotme,2)], plotme(:,sortIdx)');
        colorbar; colormap(parula(256));
        title(sprintf('%s to reward',fp_lbl)); axis square
        xlabel('time to rew (s)'); 
        ylabel('mouse'); ylim([1 nAn])
        h(y) = h.Parent; h(y).CLim = y_lims(y,:);
    %
    subplot(2, 3, 2 + ((y-1)*3)); hold on
        shadederrbar(t_rew, nanmean(plotme,2), SEM(plotme,2), clr{y});
        plot([0 0],y_lims(y,:),'--','Color',[0 0 0 0.2]);
        title(sprintf('%s to reward',fp_lbl)); axis square
        xlabel('time to rew (s)'); xticks([-1:2]);
        ylabel('%dF/F');
    %    
    subplot(2,3,3); hold on
        shadederrbar(t_rew, nanmean(a_rew_norm{y},2), SEM(a_rew_norm{y},2), clr{y}); % photometry normalized
end
subplot(2,3,3); hold on
    plot([0 0],[-1.1 1.1],'--','Color',[0 0 0 0.2]);
    title('normalized fp'); axis square
    xlabel('time to rew (s)'); xlim([-0.2 0.8]); xticks([0 0.5]);
    ylabel('fluorescence norm'); ylim([-1.1 1.1]); yticks([-1:1]);
%    
subplot(2,3,6); hold on
    plotme = lag;
    plot(plotme,'--','Color',[0 0 0 0.2]); % plot for all mice
    errorbar(nanmean(plotme,2),SEM(plotme,2),'.-k','MarkerSize',20); % plot average
    xlim([0.5 3.5]); xticks([1:3]); 
    xticklabels({'ACh peak','DA peak','ACh trough'});
    ylabel('time to reward (s)'); ylim([0 0.6]); yticks([0:0.2:0.6]);
    [~,p] = ttest(lag(:,1),lag(:,2)); % compare lag for ACh peak and DA peak
    [~,p(2)] = ttest(lag(:,2),lag(:,3)); % compare lag for DA peak and ACh trough
    title('Latency for ACh peak, DA peak, ACh trough');
    
%
movegui(gcf,'center');