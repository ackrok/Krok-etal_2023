%%
load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_figS3_dms');

corr_achda = outCorr(1).corr.corr_achda; 
shuff_achda = outCorr(1).corr.shuff_achda;
lags = outCorr(1).corr.lags;
minVal = outCorr(1).corr.minVal;
minLag = outCorr(1).corr.minLag;

f = outCorr(1).coher.f;
coher_achda = outCorr(1).coher.coher; 
coher_shuff = outCorr(1).coher.coher_shuff;
phase_achda = outCorr(1).coher.phase; 
phase_shuff = outCorr(1).coher.phase_shuff;
coherMid = outCorr(1).coher.coherMid;
phaseMid = outCorr(1).coher.phaseMid;
coherEx = outCorr(1).coher.exCoher;
phaseEx = outCorr(1).coher.exPhase;
t = outCorr(1).coher.t;

fp2phMid = outCorr(1).fpph.mid;
fp2ph = outCorr(1).fpph.fp2ph;
da2ach = outCorr(1).fpph.da2ach;
fp2phNorm = outCorr(1).fpph.fp2phNorm;

nStates = 3; nAn = 5;

% out(1).lbl = 'DMS';
% corr.corr_achda = corr_achda; corr.lags = lags; corr.shuff_achda = shuff_achda; corr.minVal = min_val; corr.minLag = min_lag;
% coher.coher = coher_achda; coher.phase = phase_achda; coher.f = f; coher.coher_shuff = coher_shuff; coher.phase_shuff = phase_shuff; coher.coherMid = coher_avg; coher.phaseMid = phase_avg;
% sig = [beh(6).FP{2}, beh(6).FP{1}];
% sig = sig - nanmean(sig);
% sig = sig(find(beh(6).time == 1200):find(beh(6).time == 1600),:);
% [ch,ph,t,f] = bz_MTCoherogram(sig(:,1),sig(:,2),'frequency',Fs,'range',[0 10],'window',10,'overlap',5,'step',5,'tapers',[3 5],'pad',0);
% coher.exCoher = ch;
% coher.exPhase = ph;
% coher.t = t;
% fpph.fp2ph = fp2ph; fpph.fp2phNorm = fp2ph_norm; fpph.da2ach = da2achph_norm; fpph.mid = mid;
% out(1).corr = corr; out(1).coher = coher; out(1).fpph = fpph;

%%
fig = figure;
fig.Position([3 4]) = [1375 1500];
clr = {'r','g','b'};

% PLOT CROSS-CORRELATION
subplot(4,3,1); hold on  
    Fs = 50;       
    for z = 1:nStates
        shadederrbar(lags/Fs, nanmean(corr_achda{z,1},2), SEM(corr_achda{z,1},2), clr{z});
    end
    shadederrbar(lags/Fs, nanmean(shuff_achda{1,2},2) - nanmean(nanmean(shuff_achda{1,2},2)),...
        nanmean(shuff_achda{1,1},2), 'k'); hold on
    plot([0 0],[-1 0.5],'--k');
    xlabel('lag from DA (s)'); xlim([-1 1]);
    ylabel('coefficient'); ylim([-0.6 0.25]); yticks([-1:0.25:1]);
    legend({'imm','mov','rew'},'Location','southeast');
    title(sprintf('DA/ACh correlation (n = %d mice)',nAn));
    axis square; set(gca,'TickDir','out');

% PLOT COHERENCE
band = [0.5 4];
subplot(4,3,5); hold on
    for z = 1:nStates % iterate over behavioral states
        shadederrbar(f, nanmean(coher_achda{z},2), SEM(coher_achda{z},2), clr{z});
    end
    plot([band(1) band(1)],[0 1],'--k'); % plot lines at 0.5, 4 Hz 
    plot([band(2) band(2)],[0 1],'--k');
    shadederrbar(f, nanmean(coher_shuff{2},2), nanmean(coher_shuff{3},2) - nanmean(coher_shuff{2},2), 'k');
    legend({'imm','loc','rew','0.5Hz','4Hz'},'Location','northeast');
    xlabel('frequency');
    ylabel('coherence'); ylim([0 1]); yticks([0:0.2:1]);
    title('coherence magnitude');
    axis square; set(gca,'TickDir','out');
    
% PLOT PHASE OFFSET
subplot(4,3,8); hold on
    for z = 1:nStates % iterate over behavioral states
        shadederrbar(f, -nanmean(rad2deg(phase_achda{z}),2), SEM(rad2deg(phase_achda{z}),2), clr{z}); 
    end
    plot([band(1) band(1)],[-180 180],'--k'); % plot lines at 0.5, 4 Hz
    plot([band(2) band(2)],[-180 180],'--k');
    shadederrbar(f, -nanmean(rad2deg(phase_shuff{2}),2), nanmean(rad2deg(phase_shuff{3}),2) - nanmean(rad2deg(phase_shuff{2}),2), 'k');
    xlabel('frequency');
    ylabel('degrees'); ylim([-180 180]); yticks([-180:90:180]);
    title('coherence phase'); 
    axis square; set(gca,'TickDir','out');

% PLOT DA to DA phase
for y = 1:2
    subplot(4,3,9+y);
    switch y; case 1; lbl = 'ACh'; case 2; lbl = 'DA'; end
    for z = 1:nStates
        shadederrbar(fp2phMid, nanmean(fp2ph{y,z},2), SEM(fp2ph{y,z},2), clr{z});
    end
    xlabel(sprintf('%s Phase',lbl)); xlim([-180 180]); xticks([-180:90:180]);
    ylabel(sprintf('%s Fluorescence (dF/F)',lbl)); ylim([-2 3])
    title(sprintf('%s fp norm to %s phase',lbl,lbl));
    axis square; set(gca,'TickDir','out');

end
% PLOT DA TO ACh phase
subplot(4,3,12);
    y = 1; % aligned to ACh phase
    sm = 5;
    for z = 1:nStates
        shadederrbar( fp2phMid, movmean(nanmean(da2ach{y,z},2),sm), movmean(SEM(da2ach{y,z},2),sm), clr{z});
        % plot( mid, movmean(nanmean(da2ach{y,z},2),sm), clr{z});
    end
    switch y; case 1; xx = 2; case 2; xx = 1; end
    plot( fp2phMid, nanmean(fp2phNorm{xx,1},2), '--k');
    switch y; case 1; lbl = 'ACh'; case 2; lbl = 'DA'; end
    xlabel(sprintf('%s Phase',lbl)); xlim([-180 180]); xticks([-180:90:180]);
    ylabel('Fluorescence (norm)'); ylim([0 1]);
    legend({'imm','mov','rew','ACh'},'Location','southeast');
    title(sprintf('fp norm to %s phase',lbl)); 
    axis square; set(gca,'TickDir','out');

jit = []; % jitter x-values for plotting
for z = 1:nStates
    j1 = z-0.25; j2 = z+0.25; jit(:,z) = j1 + (j2-j1).*rand(nAn,1); % jitter x-values for plotting
end
% PLOT STATISTICS
stats_dms = {minVal, minLag, coherMid, -rad2deg(phaseMid)};
lbl = {'coeff (pearson, r)','lag (s)','coherence','phase offset (deg)'};
lims1 = [-1 0; -0.25 0; 0 1; -180 0];
pln = [2 3 6 9];
for x = 1:length(stats_dms)
    subplot(4,3,pln(x)); hold on;
    a = stats_dms{x};
    plot(jit,a,'.k','MarkerSize',20); % plot raw values per animal
    errorbar(nanmean(a,1), SEM(a,1),'.r', 'MarkerSize', 20);
    xlim([0.5 3.5]); xticks([1:3]); xticklabels({'imm','mov','rew'}); 
    ylabel(sprintf('%s',lbl{x})); ylim(lims1(x,:));
    p = anova1(a,[],'off');
    title(sprintf('coeff: anova p = %1.2f',p));
    axis square; set(gca,'TickDir','out');
end

% PLOT EXAMPLE COHERENCE, PHASE OFFSET
sm = 0.75;
subplot(4,3,4);
a = imgaussfilt(coherEx, sm); % filter image with Gaussian smoothing kernel of X std
h = imagesc(t, f, a, [0 1]);
colorbar; colormap(jet(256));
title('Example coherence');
axis('square'); set(gca,'TickDir','out');
xlabel('Time (s)'); ylabel('Frequency (Hz)');

subplot(4,3,7);
a = imgaussfilt(rad2deg(phaseEx), sm); % filter image with Gaussian smoothing kernel of X std\
h = imagesc(t, f, a, [-180 180]);
colorbar; colormap(jet(256));
title('phase (**color axes are flipped!)');
axis('square'); set(gca,'TickDir','out');
xlabel('Time (s)'); ylabel('Frequency (Hz)');
%
movegui(gcf,'center');

%%
fig = figure;
fig.Position([3 4]) = [1375 800];

nAn = 5; nComp = 3; nStates = 3;
lbl = {'immobility','locomotion','reward'}; clr = {'r','g','b'};
lbl2 = {'DA-dls + ACh-dls','DA-dls + ACh-dms','DA-dms + ACh-dls'};
band = [0.5 4];
jit = []; % jitter x-values for plotting
for z = 1:nComp
    j1 = z-0.25; j2 = z+0.25; jit(:,z) = j1 + (j2-j1).*rand(nAn,1); % jitter x-values for plotting
end
for x = 1:3
    if ~isempty(outComp(x).coher)
    subplot(2,3,x); hold on
    for z = 1:nStates
        shadederrbar(f, nanmean(outComp(x).coher.coher{z},2), SEM(outComp(x).coher.coher{z},2), clr{z});
    end
    plot([band(1) band(1)],[0 1],'--k'); % plot lines at 0.5, 4 Hz 
    plot([band(2) band(2)],[0 1],'--k');
    shadederrbar(f, nanmean(outComp(x).coher.coher_shuff{2},2), ...
        nanmean(outComp(x).coher.coher_shuff{3},2) - nanmean(outComp(x).coher.coher_shuff{2},2), 'k');
    legend({'imm','loc','rew','0.5Hz','4Hz'});
    xlabel('Frequency (Hz)');
    ylabel('Coherence'); ylim([0 1]); yticks([0:0.2:1]);
    title(lbl2{x});
    axis square; set(gca,'TickDir','outComp');
    end
    %
    a = [outComp(1).coherMid(:,x), outComp(2).coherMid(:,x), outComp(3).coherMid(:,x)];
    subplot(2,3,x+3); hold on
    plot(jit, a, '.', 'MarkerSize', 20, 'Color', clr{x});
    errorbar([1:nComp],nanmean(a,1),SEM(a,1),'.k','MarkerSize',20);
    plot([0.5 nComp+0.5], [0.4 0.4], '--k');
    ylabel('Coherence'); ylim([0 1]);
    xlim([0.5 nComp+0.5]); xticks([1:nComp]); xticklabels(lbl2); xtickangle(45)
    title(sprintf('%s',lbl{x}));
    axis square; set(gca,'TickDir','outComp');
end
movegui(gcf,'center');

%%
% [coher_achda, phase_achda, t, f, coher_shuff, phase_shuff] = AK_coherFP(beh);
% r = [6:42]; % [~,r2] = min(abs(f - 2)); % range for 0.5-4Hz
% coher_avg = []; phase_avg = []; % initialize outCompput
% for z = 1:nStates % iterate over behavioral states
%     coher_avg(:,z) = median(coher_achda{z}(r,:)); phase_avg(:,z) = median(phase_achda{z}(r,:)); % median within frequency band
% end
% coher = struct;
% coher.coher = coher_achda; coher.phase = phase_achda; coher.f = f; 
% coher.coher_shuff = coher_shuff; coher.phase_shuff = phase_shuff; 
% outComp(x).lbl = 'DA-dms + ACh-dls';
% outComp(x).coher = coher;
% outComp(x).coherMid = coher_avg;
% outComp(x).phaseMid = phase_avg;