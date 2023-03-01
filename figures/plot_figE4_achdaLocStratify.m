load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig1_beh.mat')

%% Photometry to acceleration peak (ACh, DA)
nSplit = 3; % Divide accelerations into top, middle, bottom third
starec = cell(2,nSplit); 
for x = 1:length(beh)
    acc = getAcc(beh(x).vel); % Extract acceleration signal
    [pks,locs] = findpeaks(acc,'MinPeakProminence',0.5,'MinPeakDistance',0.5); % Location of peaks, using findpeaks function
    ev = beh(x).time(locs); % Convert peak locations to seconds
    for y = 1:length(beh(x).FP)
        fp = beh(x).FP{y}; fp = fp - nanmean(fp); Fs = beh(x).Fs;
        tmp = nan(nSplit*ceil(length(ev)/nSplit),1);
        tmp(1:length(ev)) = pks;
        [tmp,idx] = sort(tmp,'descend');
        tmp = reshape(tmp,[length(tmp)/nSplit, nSplit]);
        idx = reshape(idx,[length(idx)/nSplit, nSplit]);
        for z = 1:nSplit
            evSub = ev(idx(~isnan(tmp(:,z)),z));
            [staTmp, time] = getSTA(fp, evSub, Fs, [-1 1]);
            [staBase] = getSTA(fp, evSub, Fs, [-2 -0.5]);
            starec{y,z}(:,x) = nanmean(staTmp,2) - nanmean(staBase(:));
        end
    end
end
%%
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); nAn = length(uni);
sta = cell(2,nSplit); 
for x = 1:nAn
    ii = strcmp(rec,uni{x});
    for y = 1:size(starec,1)
        for z = 1:nSplit
            sta{y,z}(:,x) = nanmean(starec{y,z}(:,ii),2);
        end
    end
end

%%
fig = figure; fig.Position(3) = 1000;
fp_lbl = beh(1).FPnames; clr = {'g','c','k';'m','r','k'};
for y = 1:2
    subplot(1,2,y); hold on
    for z = 1:3
        shadederrbar(time, nanmean(sta{y,z},2), SEM(sta{y,z},2), clr{y,z});
    end
    plot([0 0],[-4 8],'--','Color',[0 0 0 0.5]);
    xlabel('time to peak acc (s)'); xticks([-1:0.5:1]); 
    ylabel(sprintf('%s FP (dF/F)',fp_lbl{y})); yticks([-4:2:8]);
    title(sprintf('%s stratifying acc peaks',fp_lbl{y}));
    legend({'top 1/3','middle 1/3','bottom 1/3'});
    axis square; set(gca,'TickDir','out');
end
subplot(1,2,1); ylim([-4 8])
subplot(1,2,2); ylim([-4 4])
movegui(gcf,'center');