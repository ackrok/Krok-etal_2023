%% INPUTS
nAn = length(comp(1).s);
lickWithin = 0.25;      % CHANGE, lick within this window
winRew = [-1 2];        % CHANGE, window for aligning signal to rewarded trials, in seconds
winBase = [-4 -1];
winPkDA = [100 500];    % CHANGE, window for DA peak, in milliseconds
winTrACh = [300 700];   % CHANGE, window for ACh trough, in milliseconds
xlbl = 'reward';
winAcc = [-1 1];        % CHANGE, window for aligning signal to acceleration, in seconds
NumStd = 2;             % CHANGE, for immobility trough/peak analysis
clr = {'k','b'};

%% Align photometry to reward
for z = 1:length(comp)
    beh = comp(z).s;
    rewAlign = cell(length(beh),length(beh(1).FP));
    for x = 1:length(beh)
        Fs = beh(x).Fs;
        if isempty(beh(x).reward) 
            comp(z).a_rew{x,1} = nan(151,1);
            comp(z).a_rew{x,2} = nan(151,1);
            continue; end 
        rew = beh(x).reward./Fs; % reward delivery times, in seconds
        lick = beh(x).lick./Fs; % lick times, in seconds
        [rewYes, ~, lickNew] = extractRewardedTrials(rew, lick, [0 lickWithin]);
        %rewYes: indices of rewarded trials (animal licked within specified window)
        ev = rew(rewYes); % delivery times for rewarded trials, in seconds
        for y = 1:length(beh(x).FP)
            sig = beh(x).FP{y}; % signal that will be aligned to event times
            sig = sig - nanmean(sig); % subtract mean of trace to center on zero
            [mat, time] = getSTA(sig, ev, Fs, [winRew(1), winRew(end)]);
            [matBase] = getSTA(sig, ev, Fs, [winBase(1), winBase(end)]);
            matBase = nanmean(matBase,1); % average across baseline window
            rewAlign{x,y} = mat - matBase; % save into structure
        end
    end
    rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end % mouse ID
    uni = unique(rec); % identify unique mouse ID
    comp(z).a_rew = cell(length(uni),length(beh(x).FP));
    for x = 1:length(uni)
        ii = find(strcmp(rec,uni{x}));
        for y = 1:length(beh(x).FP)
            comp(z).a_rew{x,y} = nanmean([(rewAlign{ii,y})],2); % average across multiple recordings from single mouse
        end
    end
end

%% PLOT average reward responses
if ~isfield(comp(1),'a_rew')
    error('No field a_rew - must align photometry to reward delivery before proceeding.');
end
fig = figure; fig.Position(3) = 1000;
for z = 1:2
    a_rew = comp(z).a_rew; % extract reward-aligned data from structure
    for y = 1:size(a_rew,2)
        a_mat = []; % initialize temporary output matrix
        for x = 1:size(a_rew,1)
            a_mat(:,x) = nanmean(a_rew{x,y},2); % extract reward-aligned data from structure and average across all rewards
        end
        subplot(1,2,y); hold on
        plot(time, a_mat, clr{z});
        % shadederrbar(time, nanmean(a_mat,2), SEM(a_mat,2), clr{z}); % plot average across trials
        ylabel(sprintf('%s amplitude (dF/F)',comp(1).s(1).FPnames{y}));
        xlabel(sprintf('time to %s (s)',xlbl));
        title(sprintf('%s (n = %d): %s vs %s',...
            comp(1).s(1).FPnames{y},...
            size(a_mat(:,~isnan(a_mat(1,:))),2),...
            comp(1).lbl,comp(2).lbl)); 
        axis square
    end
end
movegui(gcf,'center');

%% PLOT example reward responses
x = 4; % example
if ~isfield(comp(1),'a_rew')
    error('No field a_rew - must align photometry to reward delivery before proceeding.');
end
fig = figure; fig.Position(3) = 1000;
for z = 1:2
    a_rew = comp(choice(z)).a_rew; % extract reward-aligned data from structure
    for y = 1:size(a_rew,2)
        subplot(1,2,y);
        shadederrbar(time, nanmean(a_rew{x,y},2), SEM(a_rew{x,y},2), clr{choice(z)}); % plot average across trials
        ylabel(sprintf('%s amplitude (dF/F)',comp(1).s(1).FPnames{y}));
        xlabel(sprintf('time to %s (s)',xlbl));
        title(sprintf('%s (example): %s vs %s',...
            comp(1).s(1).FPnames{y},...
            comp(choice(1)).inf,comp(choice(2)).inf)); 
        axis square
    end
end
movegui(gcf,'center');

%% Find AMPLITUDE and LATENCY of reward response
if ~isfield(comp(1),'a_rew')
    error('No field a_rew - must align photometry to reward delivery before proceeding.');
end
fig = figure; fig.Position(3) = 1000;
lag = {}; val = {}; % initialize output
for z = 1:2
    a_rew = comp(z).a_rew; % extract reward-aligned data from structure
    Fs = comp(z).s(1).Fs; % sampling frequency
    y = 1; % analyzing ACh reward response
    win = winTrACh; % range for ACh trough, in milliseconds
    r = (win/(1000/Fs) + find(time == 0));
    for x = 1:size(a_rew,1)
        pull = a_rew{x,y}; % extract matrix of reward-aligned signals during post-infusion window into workspace
        pull = pull(r(1):r(2),:); % extract only segment of time specified as window for reward-related trough
        [a,b] = min(pull); % find local MINIMUM within range
        c = time(b + r(1) - 1); % convert index to seconds
        c(isnan(a)) = nan; 
        val{y}(x,z) = nanmean(a); % save amplitude of minima
        lag{y}(x,z) = nanmean(c).*1000; % save lag at minima, in milliseconds
    end
    val{y}(val{y} == 0) = nan; lag{y}(lag{y} == 0) = nan;
    
    y = 2; % analyzing DA reward response
    win = winPkDA; % range for ACh trough, in milliseconds
    r = (win/(1000/Fs) + find(time == 0)); 
    for x = 1:size(a_rew,1)
        pull = a_rew{x,y}; % extract matrix of reward-aligned signals during post-infusion window into workspace
        pull = pull(r(1):r(2),:); % extract only segment of time specified as window for reward-related trough
        [a,b] = max(pull); % find local MAXIMUM within range
        c = time(b + r(1) - 1); % convert index to seconds
        c(isnan(a)) = nan; 
        val{y}(x,z) = nanmean(a); % save amplitude of minima
        lag{y}(x,z) = nanmean(c).*1000; % save lag at minima, in milliseconds
    end
    val{y}(val{y} == 0) = nan; lag{y}(lag{y} == 0) = nan;
    
    for y = 1:2
        subplot(1,2,y); hold on
        plot(lag{y}(:,z), val{y}(:,z), '.', 'MarkerSize', 20, 'Color', clr{z}); % plot individual data points for each mouse
        errorbar(nanmean(lag{y}(:,z)), nanmean(val{y}(:,z)), ... % plot average across all mice
        SEM(val{y}(:,z),1), SEM(val{y}(:,z),1), ... % error bars vertical for amplitude
        SEM(lag{y}(:,z),1), SEM(lag{y}(:,z),1), ... % error bars horizontal for latency
        '.', 'MarkerSize', 20, 'Color', clr{z});
        xlabel(sprintf('time to %s (s)',xlbl)); xlim([win(1) win(2)]);
    end
end
y = 1; subplot(1,2,y); % ACh reward response subplot
    xlabel('time to reward (s)'); xlim([winTrACh(1) winTrACh(2)]);
    ylabel('ACh trough amplitude (%dF/F)'); ylim([-4 0]); yticks([-8:1:0]);
    [~,p] = ttest2(lag{y}(:,1),lag{y}(:,2)); % statistical test: paired t-test
    [~,p(2)] = ttest2(val{y}(:,1),val{y}(:,2)); % statistical test: paired t-test
    title(sprintf('%s vs %s (lag %1.2f, val %1.2f)',comp(1).lbl,comp(2).lbl,p(1),p(2)));
    axis square
y = 2; subplot(1,2,y); % DA reward response subplot
    xlabel('time to reward (s)'); xlim([winPkDA(1) winPkDA(2)]);
    ylabel('DA peak amplitude (%dF/F)'); ylim([0 10]); yticks([0:1:10]);
    [~,p] = ttest2(lag{y}(:,1),lag{y}(:,2)); % statistical test: paired t-test
    [~,p(2)] = ttest2(val{y}(:,1),val{y}(:,2)); % statistical test: paired t-test
    title(sprintf('%s vs %s (lag %1.2f, val %1.2f)',comp(1).lbl,comp(2).lbl,p(1),p(2)));
    axis square
movegui(gcf,'center');