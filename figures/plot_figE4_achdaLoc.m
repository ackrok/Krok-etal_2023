load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig1_beh.mat')

%% Extract standard deviation (ACh, DA)
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
sig_an = cell(1,2);
for x = 1:nAn
    ii = find(strcmp(rec,uni{x}));
    for y = 1:2
        sig_an{y}(x,:) = nanmean(sig{y}(ii,:),1); % average across multiple recordings from same mouse
    end
end
fprintf('Extracting standard deviation DONE\n');

fig = figure; fig.Position(3) = 1000;
clr = {'g','m'}; lbl = beh(1).FPnames;
for y = 1:2
    subplot(1,2,y); hold on
    plot(sig_an{y}', '.-', 'Color', [0 0 0 0.1], 'MarkerSize', 20);
    errorbar([0.75 2.25], nanmean(sig_an{y}), SEM(sig_an{y},1), '.', 'MarkerSize', 20, 'Color', clr{y});
    xticks([1 2]); xticklabels({'imm','loc'}); 
    ylabel('standard deviation'); ylim([0 6]);
    [~,p] = ttest(sig_an{y}(:,1),sig_an{y}(:,2));
    title(sprintf('%s STD (p = %1.4f) n = %d',lbl{y},p,nAn)); axis square
end
movegui(gcf,'center');

%% Acceleration to ACh/DA phase (locomotion, immobility)
fig = figure; fig.Position(3) = 1000;

behAcc = beh;
for x = 1:length(behAcc); behAcc(x).FP{2} = getAcc(behAcc(x).vel); behAcc(x).FPnames{2} = 'acceleration'; end
[~, ~, ~, mid, acc2achph] = AK_fp2phase(behAcc);
subplot(1,2,1); hold on
y = 1; clr = {'k','g'};
for z = 1:2; shadederrbar(mid, nanmean(acc2achph{y,z},2), SEM(acc2achph{y,z},2), clr{z}); end
xlabel('ACh phase'); xlim([-180 180]); 
ylabel('Acceleration (cm/s^2)'); yticks([-0.1:0.1:0.2]);
title('Locomotion/immobility to ACh phase');
axis square; set(gca,'TickDir','out');

behAcc = beh;
for x = 1:length(behAcc)
    behAcc(x).FP{1} = behAcc(x).FP{2}; behAcc(x).FPnames{1} = behAcc(x).FPnames{2};  % swap DA and ACh photometry 
    behAcc(x).FP{2} = getAcc(behAcc(x).vel); behAcc(x).FPnames{2} = 'acceleration'; % {1} DA, {2} acceleration
end
[~, ~, ~, mid, acc2achph] = AK_fp2phase(behAcc);
subplot(1,2,2); hold on
y = 1; clr = {'k','m'};
for z = 1:2; shadederrbar(mid, nanmean(acc2achph{y,z},2), SEM(acc2achph{y,z},2), clr{z}); end
xlabel('DA phase'); xlim([-180 180]); xticks([-180:90:180]);
ylabel('Acceleration (cm/s^2)'); yticks([-0.1:0.1:0.2]);
title('Locomotion/immobility to DA phase');
axis square; set(gca,'TickDir','out');
movegui(gcf,'center');

%% Photometry to acceleration peak (ACh, DA)
starec = cell(1,2); 
for x = 1:length(beh)
    acc = getAcc(beh(x).vel); % Extract acceleration signal
    [~,locs] = findpeaks(acc,'MinPeakProminence',0.5,'MinPeakDistance',0.5); % Location of peaks, using findpeaks function
    ev = beh(x).time(locs); % Convert peak locations to seconds
    for y = 1:length(beh(x).FP)
        fp = beh(x).FP{y}; fp = fp - nanmean(fp); Fs = beh(x).Fs;
        [staTmp, time] = getSTA(fp, ev, Fs, [-1 1]);
        [staBase] = getSTA(fp, ev, Fs, [-2 -0.5]);
        starec{y}(:,x) = nanmean(staTmp,2) - nanmean(staBase(:));
    end
end
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); nAn = length(uni);
sta = cell(1,2); 
for x = 1:nAn
    ii = strcmp(rec,uni{x});
    for y = 1:length(starec)
        sta{y}(:,x) = nanmean(starec{y}(:,ii),2);
    end
end
m = max(sta{1}); [~,sortIdx] = sort(m); % sort by maximal ACh to acc

%
fig = figure; fig.Position(3) = 1375;
fp_lbl = beh(1).FPnames; clr = {'g','m'};
y_lims = [-2 10; -2 4];
for y = 1:length(sta)
    subplot(1,3,y); hold on
    h = imagesc(time, [1:nAn], sta{y}(:,sortIdx)');
    colorbar; colormap(parula(256));
    title(sprintf('%s to peak acc.',fp_lbl{y})); axis square
    xlabel('time to peak acc. (s)'); xticks([-1:0.5:1]);
    ylabel('mouse'); ylim([0 nAn+1])
    h(y) = h.Parent; h(y).CLim = y_lims(y,:);
    
    subplot(1,3,3); hold on
    shadederrbar(time, nanmean(sta{y},2), SEM(sta{y},2), clr{y});
    plot([0 0],[-2 5],'--','Color',[0 0 0 0.2]);
    xlabel('time to peak acc. (s)'); xticks([-1:0.5:1]); 
    ylabel('FP (% dF/F)');
    axis square; set(gca,'TickDir','out');
end

movegui(gcf,'center');