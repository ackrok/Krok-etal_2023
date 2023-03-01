load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig1_beh');
load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_figS7_loc');

%% ACh photometry: locomotion/immobility standard deviation
idxStates = extractBehavioralStates(beh); % identify indices pertaining to different behavioral states
sig = cell(1,2);
for x = 1:length(beh)
    fp = [beh(x).FP{1}, beh(x).FP{2}]; % extract photometry signal
    fp = fp - nanmean(fp);
    for y = 1:2
        sig{y}(x,1) = nanstd(fp(idxStates{x,1},y)); % standard deviation of photometry signal during immobility
        sig{y}(x,2) = nanstd(fp(idxStates{x,2},y)); % standard deviation of photometry signal during locomotion
    end
end
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); nAn = length(uni); % number of unique animals
stdFP = cell(1,2);
for x = 1:nAn
    ii = find(strcmp(rec,uni{x}));
    for y = 1:2
        stdFP{y}(x,:) = nanmean(sig{y}(ii,:),1); % average across multiple recordings from same mouse
    end
end

%% ACh photometry: locomotion/immobility acceleration to ACh phase
behAcc = beh;
for x = 1:length(behAcc); behAcc(x).FP{2} = getAcc(behAcc(x).vel); behAcc(x).FPnames{2} = 'acceleration'; end
[~, ~, ~, mid, acc2achph] = AK_fp2phase(behAcc);

%% ACh photometry: locomotion/immobility ratio
idxStates = extractBehavioralStates(beh); % identify indices pertaining to different behavioral states
delta = nan(length(beh),1); % initialize output vector
for x = 1:length(beh)
    fp = beh(x).FP{1};
    tmp = [nanmean(fp(idxStates{x,1})), nanmean(fp(idxStates{x,2}))]; % average photometry signal during immobility, locomotion
    delta(x) = (tmp(2) - tmp(1))/tmp(1); % delta change from immobility to locomotion
end
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); % unique mouse IDs
deltaFP = nan(length(uni),1);
for x = 1:length(uni) % iterate over all unique mouse IDs
    ii = find(strcmp(rec, uni{x})); % identify recordings for one mouse
    deltaFP(x) = nanmean(delta(ii));
end

%% pCIN firing rate: locomotion/immobility ratio
delta = nan(length(units),1); % initialize output vector 
for x = 1:length(units)
    b = find(strcmp({behUnits.rec},units(x).rec)); 
    if isempty(behUnits(b).on); continue; end
    st = units(x).st;
    tmp = 1/mean(diff(extractEventST(st,behUnits(b).onRest,behUnits(b).offRest,0))); % firing rate during immobility
    tmp(2) = 1/mean(diff(extractEventST(st,behUnits(b).on,behUnits(b).off,0))); % firing rate during locomotion
    delta(x) = (tmp(2) - tmp(1))/tmp(1);
end
rec = {}; for x = 1:length(units); rec{x} = strtok(units(x).rec,'_'); end
uni = unique(rec); % unique mouse IDs
deltaFR = nan(length(uni),1);
for x = 1:length(uni) % iterate over all unique mouse IDs
    ii = find(strcmp(rec, uni{x})); % identify recordings for one mouse
    deltaFR(x) = nanmean(delta(ii));
end

%% PLOT
fig = figure; fig.Position(3) = 1375;
subplot(1,3,1); hold on
    y = 1;
    plot(stdFP{y}', '.-', 'Color', [0 0 0 0.1], 'MarkerSize', 20);
    errorbar([0.75 2.25], nanmean(stdFP{y}), SEM(stdFP{y},1), '.g', 'MarkerSize', 20);
    xticks([1 2]); xticklabels({'immobility','locomotion'}); 
    ylabel('Standard deviation'); ylim([0 6]);
    [~,p] = ttest(stdFP{y}(:,1),stdFP{y}(:,2));
    title(sprintf('STD ACh: p = %1.4f',p));
    axis square; set(gca,'TickDir','out');
subplot(1,3,2); hold on
    y = 1; clr = {'k','g'};
    for z = 1:2; shadederrbar(mid, nanmean(acc2achph{y,z},2), SEM(acc2achph{y,z},2), clr{z}); end
    xlabel('ACh phase'); xlim([-180 180]); 
    ylabel('Acceleration (cm/s^2)'); yticks([-0.1:0.1:0.2]);
    title('Locomotion/immobility to ACh phase');
    axis square; set(gca,'TickDir','out');
subplot(1,3,3); hold on
    plot([0.5 2.5],[1 1],'--k');
    nAn = length(deltaFP);
    j1 = 0.8; j2 = 1.2; jit = j1 + (j2-j1).*rand(nAn,1); % jitter
    plot(jit, 1+deltaFP, '.g', 'MarkerSize', 20); % Plot area under curve power data points for each recording
    errorbar(1, 1+nanmean(deltaFP), SEM(deltaFP,1), '.k', 'MarkerSize', 20); % Plot area under curve power averaged
    nAn = length(deltaFR);
    j1 = 1.8; j2 = 2.2; jit = j1 + (j2-j1).*rand(nAn,1); % jitter
    plot(jit, 1+deltaFR, '.b', 'MarkerSize', 20); % Plot area under curve power data points for each recording
    errorbar(2, 1+nanmean(deltaFR), SEM(deltaFR,1), '.k', 'MarkerSize', 20); % Plot area under curve power averaged
    xlim([0.5 2.5]); xticks([1 2]); xticklabels({'ACh','pCIN'});
    ylabel('Locomotion/immobility ratio'); yticks([0.5:0.5:2]);
    [~,p] = ttest2(deltaFP,deltaFR);
    title(sprintf('ACh vs CIN: p = %1.4f',p));
    axis('square'); set(gca,'TickDir','out');
movegui(gcf,'center');

%% PLOT AVERAGE
fig = figure; fig.Position([3 4]) = [1000 800];
lims1 = [-1 2; -8 14]; lims2 = [-0.25 0.5; -3 7]; lims3 = [-0.5:0.25:0.5; -2:2:6];
clr = {'b','g'};
for x = 1:length(out)
    subplot(2,2,x); hold on
        [~, ii] = sort(max(out(x).delta)); % sort in ascending order 
        ii = fliplr(ii);
        h = imagesc(out(x).time, [1:out(x).n], out(x).delta(:,ii)', [0.6 2.2]);
        colorbar; colormap(jet(256));
        xlabel('Time to peak acc. (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
        ylabel('unit/mouse');
        h(x) = h.Parent; h(x).CLim = lims1(x,:);
        title(sprintf('%s',out(x).lbl));
        axis square; set(gca,'TickDir','out');

    subplot(2,2,x+length(out)); hold on
        sm = 5;
        plot([-1 1],[0 0],'--k');
        plot([0 0],[lims2(x,1) lims2(x,2)],'--k');
        shadederrbar(out(x).time, nanmean(movmean(out(x).delta,1,1),2), SEM(movmean(out(x).delta,1,1),2), clr{x}); 
        xlabel('Time to peak acc. (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
        ylabel(sprintf('%s',out(x).ylbl)); ylim(lims2(x,:)); yticks(lims3(x,:));
        title(sprintf('%s (n = %d)',out(x).lbl,out(x).n));
        axis('square'); set(gca,'TickDir','out');
end
movegui(gcf,'center');