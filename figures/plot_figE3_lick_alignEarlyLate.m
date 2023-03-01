load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig1_reward_beh.mat')

%%
seg = 4; % how many segments to divide recording into
lickWithin = 1; % lick within this window
winRew = [-1 2]; % window for aligning signal to rewarded trials
winBase = [-4 -1]; % window for baseline
staRec = cell(3,2); % initialize output

for x = 1:length(beh)
    Fs = beh(x).Fs;
    rew = beh(x).reward(:)./Fs; % reward delivery times, in seconds
    lick = beh(x).lick(:)./Fs; % lick times, in seconds
    [rewYes, rewNo, lickNew] = extractRewardedTrials(rew, lick, [0 lickWithin]);
    
    %% early vs late
    segment = floor(length(rew)/seg); % How many rewards included in each segment
    rewEarly = rewYes(rewYes < segment); % rewarded trials during first segment of reward deliveries
    rewLate = rewYes(rewYes > (length(rew)-segment)); % rewarded trials during last segment of reward deliveries

    %%
    ev = cell(3,1);
    ev{1} = rew(rewEarly); ev{2} = rew(rewLate);
    ev{3} = randi(length(beh(x).FP{1}),length(rew),1);
    
    %% aligning photometry
    for y = 1:length(beh(x).FP)
        fp = beh(x).FP{y}; fp = fp - nanmean(fp);
        for z = 1:length(ev)
            [mat, time] = getSTA(fp, ev{z}, Fs, [winRew(1), winRew(end)]);
            [matBase] = getSTA(fp, ev{z}, Fs, [winBase(1), winBase(end)]);
            matBase = nanmean(matBase,1); % Average across baseline window
            staRec{z,y}(:,x) = nanmean(mat - matBase, 2);
        end
    end
end

rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); nAn = length(uni); % Number of unique animals
staAn = cell(3,2);
for x = 1:nAn
    ii = find(strcmp(rec,uni{x}));
    for y = 1:2
        for z = 1:3
            staAn{z,y}(:,x) = nanmean(staRec{z,y}(:,ii),2); % average across multiple recordings from same mouse
        end
    end
end
fprintf('DONE\n');

%% Peak amplitude and latency for all mice, early vs late
nAn = 14;
winPkDA = [100 500];    % CHANGE, window for DA peak, in milliseconds
winTrACh = [300 700];   % CHANGE, window for ACh trough, in milliseconds
lag = {}; val = {}; % initialize output
for z = 1:2
    Fs = beh(x).Fs; % sampling frequency
    y = 1; % analyzing ACh reward response
    win = winTrACh; % range for ACh trough, in milliseconds
    r = (win/(1000/Fs) + find(time == 0));
    pull = staRec{z,y}; % extract matrix of reward-aligned signals during post-infusion window into workspace
    pull = pull(r(1):r(2),:); % extract only segment of time specified as window for reward-related trough
    [a,b] = min(pull); % find local MINIMUM within range
    c = time(b + r(1) - 1); % convert index to seconds
    c(isnan(a)) = nan; 
    val{y}(:,z) = a; % save amplitude of minima
    lag{y}(:,z) = c.*1000; % save lag at minima, in milliseconds
    val{y}(val{y} == 0) = nan; lag{y}(lag{y} == 0) = nan;
    
    y = 2; % analyzing DA reward response
    win = winPkDA; % range for ACh trough, in milliseconds
    r = (win/(1000/Fs) + find(time == 0)); 
    pull = staRec{z,y}; % extract matrix of reward-aligned signals during post-infusion window into workspace
    pull = pull(r(1):r(2),:); % extract only segment of time specified as window for reward-related trough
    [a,b] = max(pull); % find local MINIMUM within range
    c = time(b + r(1) - 1); % convert index to seconds
    c(isnan(a)) = nan; 
    val{y}(:,z) = a; % save amplitude of minima
    lag{y}(:,z) = c.*1000; % save lag at minima, in milliseconds
    val{y}(val{y} == 0) = nan; lag{y}(lag{y} == 0) = nan;
end

%% Plot photometry to rewarded trials
fig = figure; fig.Position(3) = 1375;
clr = {'g','m';'b','r'};
fig_lbl = {'rewarded trials (early vs late)'};

subplot(1,3,1); hold on;
    for z = 1:2
        for y = 1:2
            shadederrbar(time, nanmean(staRec{z,y},2), SEM(staRec{z,y},2), clr{z,y}); % plot average across trials
        end
    end
    y = 1; mm = abs(nanmean(SEM(staRec{3,y},2))) + abs(nanmean(nanmean(staRec{3,y},2)));
    shadederrbar(time, zeros(length(time),1), mm.*ones(length(time),1), 'k'); 
    plot([0 0],[-5 10],'-','Color',[0 0 0 0.8]);
    xlabel('time to reward (s)'); xticks([-1:0.5:2]);
    ylabel('FP (%dF/F)'); yticks([-4:4:12]);
    legend({'earlyACh','earlyDA','lateACh','lateDA'});
    title(sprintf('%s ',fig_lbl{1}));
    axis square; set(gca,'TickDir','out');

%%
winTrPk = [winTrACh; winPkDA]; clr = {'k','r'}; lbl_fp = {'ACh','DA'};
for y = 1:2
    subplot(1,3,y+1); hold on
    for z = 1:2
        for x = 1:size(val{y},1)
            plot(lag{y}(x,:), val{y}(x,:), '-', 'Color', [0 0 0 0.1],'HandleVisibility','off');
        end
        plot(lag{y}(:,z), val{y}(:,z), '.', 'MarkerSize', 20, 'Color', clr{z}); % plot individual data points for each mouse
        errorbar(nanmean(lag{y}(:,z)), nanmean(val{y}(:,z)), ... % plot average across all mice
        SEM(val{y}(:,z),1), SEM(val{y}(:,z),1), ... % error bars vertical for amplitude
        SEM(lag{y}(:,z),1), SEM(lag{y}(:,z),1), ... % error bars horizontal for latency
        '.', 'MarkerSize', 20, 'Color', clr{z},'HandleVisibility','off');
    end
    xlabel('time to reward (s)'); xlim([winTrPk(y,:)]);
    ylabel('amplitude (%dF/F)'); % ylim([0 15]);
    legend({'early','late'},'Location','northeast');
    [~,p] = ttest(lag{y}(:,1),lag{y}(:,2)); % statistical test: paired t-test
    [~,p(2)] = ttest(val{y}(:,1),val{y}(:,2)); % statistical test: paired t-test
    title(sprintf('%s early vs late (lag %1.2f, val %1.2f)',lbl_fp{y},p(1),p(2)));
    axis square
end
movegui(gcf,'center');

%% UNITY LINE
winTrPk = [winTrACh; winPkDA]; clr = {'k','r'}; lbl_fp = {'ACh','DA'};
for y = 1:2
    subplot(1,3,y+1); hold on
    plot(val{y}(:,1), val{y}(:,2), '.k', 'MarkerSize', 20);
    switch y; case 1; plot([-10:0],[-10:0],'-r');
        case 2; plot([0:20],[0:20],'-r'); end
    xlabel(sprintf('%s early amplitude',lbl_fp{y})); 
    ylabel(sprintf('%s late amplitude',lbl_fp{y})); 
    [~,p] = ttest(lag{y}(:,1),lag{y}(:,2)); % statistical test: paired t-test
    [~,p(2)] = ttest(val{y}(:,1),val{y}(:,2)); % statistical test: paired t-test
    title(sprintf('%s early vs late (lag %1.2f, val %1.2f)',lbl_fp{y},p(1),p(2)));
    axis square; set(gca,'TickDir','out');
end
movegui(gcf,'center');