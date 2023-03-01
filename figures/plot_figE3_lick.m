load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig1_reward_beh.mat')

%%
lickWithin = 1; % lick within this window
winRew = [-1 2]; % window for aligning signal to rewarded trials
winBase = [-4 -1]; % window for baseline
staRec = cell(5,2); % initialize output

for x = 1:length(beh)
    Fs = beh(x).Fs;
    rew = beh(x).reward(:)./Fs; % reward delivery times, in seconds
    lick = beh(x).lick(:)./Fs; % lick times, in seconds
    [rewYes, rewNo, lickNew] = extractRewardedTrials(rew, lick, [0 lickWithin]);
    
    %% 1st rewarded lick
    bin = 1/1000; 
    peth = getClusterPETH(lickNew, rew(rewYes), bin, [0 1]); % Lick aligned to rewarded trials in 1 ms bins
    cts = peth.cts{1}; % Lick counts in 1ms bins for each reward trial
    [~, lickFirst] = max(cts~=0, [], 1); % Find first non-zero index for each trial
    lickFirst = lickFirst(:).*bin + rew(rewYes); % First lick after delivery for rewarded trials
    
    %% all non-rewarded licks
    lickRew = extractEventST(lickNew, rew + winRew(1), rew + winRew(2), 1); % Odentify licks during reward period
    lickNonRew = lickNew(~ismember(lickNew, lickRew)); % Exclude rewarded licks
    
    %%
    ev = cell(5,1);
    ev{1} = rew(rewYes); ev{2} = rew(rewNo); 
    ev{3} = lickFirst; ev{4} = lickNonRew;
    ev{5} = randi(length(beh(x).FP{1}),length(rew),1);
    
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
staAn = cell(4,2);
for x = 1:nAn
    ii = find(strcmp(rec,uni{x}));
    for y = 1:2
        for z = 1:4
            staAn{z,y}(:,x) = nanmean(staRec{z,y}(:,ii),2); % average across multiple recordings from same mouse
        end
    end
end
fprintf('DONE\n');
        
%% Plot photometry to rewarded trials
fig = figure; fig.Position([3 4]) = [1000 800];
clr = {'g','m'};
fig_lbl = {'rewarded trials (lick within 1s)','non-rewarded trials (no lick within 1s)','1st rewarded lick','non-rewarded licks'};
x_lbl = {'reward','reward','lick','lick'};
for z = 1:4
    sp(z) = subplot(2,2,z); hold on
    for y = 1:2
        shadederrbar(time, nanmean(staRec{z,y},2), SEM(staRec{z,y},2), clr{y}); % plot average across trials
    end
    y = 1; mm = abs(nanmean(SEM(staRec{5,y},2))) + abs(nanmean(nanmean(staRec{5,y},2)));
    shadederrbar(time, zeros(length(time),1), mm.*ones(length(time),1), 'k'); 
    plot([0 0],[-5 10],'-','Color',[0 0 0 0.8]);
    xlabel(sprintf('time to %s (s)',x_lbl{z})); xticks([-1:0.5:2]);
    ylabel('FP (%dF/F)');
    title(sprintf('%s ',fig_lbl{z}));
    axis square; set(gca,'TickDir','out');
end
linkaxes(sp,'y');
movegui(gcf,'center');
    