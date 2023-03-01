load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig1_reward_beh.mat')

%% Extract data and process lick and reward trials
out = struct;
window = [0 1]; % Window for reward to be collected
for x = 1:length(beh)
    out(x).rec = beh(x).rec;
    Fs = beh(x).Fs;
    rew = beh(x).reward(:)./Fs; % reward delivery times, in seconds
    lick = beh(x).lick(:)./Fs; % lick times, in seconds
    [rewYes, rewNo, lickNew] = extractRewardedTrials(rew, lick, window);
    
%% Extract licks, correct for timing
    lick = lickNew; % Overwrite lick vectorr
    out(x).lick_new = lick; % Corrected licks, in seconds
    
%% Identify rewarded trials
    out(x).delivery = rew; % Reward delivery times, in seconds
    ev = rew(rewYes); % Event time is rewarded trials
    out(x).rew_no = rewNo; % Index of deliveries where animal did not lick to receive reward
    out(x).rew_yes = rewYes; % Index of delivieries where animal DID lick

%% Onset of lick bouts for each rewarded trial
    bin = 1/1000;
    peth = getClusterPETH(lick, ev, bin, window); % PETH: lick aligned to rewarded trials in 1 ms bins
    cts = peth.cts{1}; % Lick counts in 1ms bins for each reward trial
    [~, lick_first] = max(cts~=0, [], 1); % Find first non-zero index for each trial
%     [~, lick_last] = max(flipud(cts)~=0, [], 1); % Find last non-zero index
%     lick_last = window(2)/bin - (lick_last - 1); % Flip back because last lick was determined from flipped matrix
    out(x).rew_lick = lick_first(:).*bin + ev(:); % First lick after delivery for rewarded trials
    tmp = nan(length(out(x).delivery),1); % Initialize output vector
    tmp(out(x).rew_yes) = out(x).rew_lick;
    out(x).rew_lick = tmp; % Include NaNs for non-rewarded trials
    
%% Number of licks per reward period (0:2s)
    peth = getClusterPETH(lick, ev, 2, [0 2]); % PETH: lick aligned to rewarded trials
    cts = peth.cts{1}(:); % Lick counts for each reward trial across 2 second period post-delivery
    out(x).lick_num = nan(length(out(x).delivery),1);
    out(x).lick_num(out(x).rew_yes) = cts;
    
end
fprintf('Done re-processing lick and reward trials. Window: %d to %1.2f seconds \n',window(1),window(2));

%% PLOT - Timing of non-rewarded trials
rewNoBin = [];
bin = 0.05; edges = [0:bin:1]; mid = edges(2:end)-bin/2;
for x = 1:length(out)
    b = (out(x).rew_no)./length(out(x).delivery); % Proportion of the way through each trial that non-rewarded trial occurs (number of non-rewarded trials divided by the total number of trials)
    n = histcounts(b, edges, 'Normalization', 'probability');
    rewNoBin(x,:) = n;
end
rewNoBin(:,1) = nan;

figure; hold on
bar(mid, nanmean(rewNoBin,1));
errorbar(mid, nanmean(rewNoBin,1), SEM(rewNoBin,1), SEM(rewNoBin,1), 'k', 'LineStyle', 'none');
xlabel('progression through trial'); xticks([0:0.25:1]);
ylabel('probability of non-rewarded trial'); yticks([0:0.1:0.4]);
a = rewNoBin(:,[2:4]); b = rewNoBin(:,[18:20]);
[~,p] = ttest(a(:),b(:));
title(sprintf('Timing of non-rewarded trials (lick within %1.2f s)? \n bin2-4 vs bin18-20: p = %1.3f',window(2),p));

%% Median delivery to lick latency (subplot 1/3)
seg = 4; % How many segments to divide recording into
lick_num = [];
for x = 1:length(out)
    segment = floor(length(out(x).delivery)/seg); % How many rewards included in each segment
    lick0 = out(x).rew_lick - out(x).delivery; % Delivery-to-lick latency
    % lick0(isnan(lick0)) = window(2); % replace NaN's with maximum
    early = lick0(1:segment); % EARLY reward delivery
    late = lick0(end-segment+1:end); % LATE reward delivery
    lick_num(x,:) = [median(early,'omitnan'), median(late,'omitnan')]; % MEDIAN lick timing
end

fig = figure; fig.Position(3) = 1375;
subplot(1,3,1); hold on
pull = lick_num.*1000; 
lbl2 = {sprintf('early delivery (1/%d)',seg),sprintf('late delivery (%d/%d)',seg,seg)};
plot([1;2].*ones(2, length(pull)), pull', '.-', 'MarkerSize', 20, 'Color', [0 0 0 0.2]);
errorbar([0.75 2.25], nanmean(pull), SEM(pull,1), '.', 'MarkerSize', 20, 'Color', 'r');
xlim([0.5 2.5]); xticks([1 2]); xticklabels(lbl2); 
ylabel('median delivery-to-lick latency (ms)'); ylim([0 350]); yticks([0:100:1000]);
[~,p] = ttest(pull(:,1),pull(:,2));
title(sprintf('delivery-to-lick (win %d-%1.2f) (p = %1.3f)',window(1),window(2),p)); 
axis square; set(gca,'TickDir','out');
% fprintf('Delivery-to-lick (median) latency: %d +/- %d ms (early) vs %d +/- %d ms (late) \n   for "rewarded" trial -- lick within %d ms \n',...
%     round(nanmean(pull(:,1))), round(SEM(pull(:,1),1)), round(nanmean(pull(:,2))), round(SEM(pull(:,2),1)), window(2)*1000);
movegui(gcf,'center');

%% Probability of rewarded trial (subplot 2/3)
a = rewNoBin(:,[2:5]); a = nanmean(a,2); a = 1-a;
b = rewNoBin(:,[17:20]); b = nanmean(b,2); b = 1-b;

subplot(1,3,2); hold on
pull = [a,b].*100; 
plot([1;2].*ones(2, length(pull)), pull', '.-', 'MarkerSize', 20, 'Color', [0 0 0 0.2]);
errorbar([0.75 2.25], nanmean(pull), SEM(pull,1), '.', 'MarkerSize', 20, 'Color', 'r');
xlim([0.5 2.5]); xticks([1 2]); xticklabels(lbl2); 
ylabel('probability of rewarded trial'); ylim([0 100]); yticks([0:25:100]);
[~,p] = ttest(pull(:,1),pull(:,2));
title(sprintf('probability of rewarded trial (p = %1.3f)',p));
axis square; set(gca,'TickDir','out');

%% Number of licks per reward (subplot 3/3)
seg = 4; % How many segments to divide recording into
lickNum = [];
for x = 1:length(out)
    segment = floor(length(out(x).delivery)/seg); % How many rewards included in each segment
    early = out(x).lick_num(1:segment); % EARLY reward delivery
    late = out(x).lick_num(end-segment+1:end); % LATE reward delivery
    lickNum(x,:) = [median(early,'omitnan'), median(late,'omitnan')]; % MEDIAN lick timing
end

subplot(1,3,3); hold on
pull = lickNum; 
plot([1;2].*ones(2, length(pull)), pull', '.-', 'MarkerSize', 20, 'Color', [0 0 0 0.2]);
errorbar([0.75 2.25], nanmean(pull), SEM(pull,1), '.', 'MarkerSize', 20, 'Color', 'r');
xlim([0.5 2.5]); xticks([1 2]); xticklabels(lbl2); 
ylabel('Number of licks per reward (in 2s)'); ylim([0 20]); yticks([0:5:20]);
[~,p] = ttest(pull(:,1),pull(:,2));
title(sprintf('Number of licks per reward (p = %1.3f)',p));
axis square; set(gca,'TickDir','out');
