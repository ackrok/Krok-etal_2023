%% INPUTS
load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig3_corr');

corr_achda = out(1).corr.corr_achda; 
shuff_achda = out(1).corr.shuff_achda;
lags = out(1).corr.lags;
minVal = out(1).corr.minVal./nanmean(out(2).corr.minVal);
minLag = out(1).corr.minLag;

f = out(1).coher.f;
coher_achda = out(1).coher.coher; coher_shuff = out(1).coher.coher_shuff;
phase_achda = out(1).coher.phase; phase_shuff = out(1).coher.phase_shuff;
coherMid = out(1).coher.coherMid./nanmean(out(2).coher.coherMid);
phaseMid = out(1).coher.phaseMid./nanmean(out(2).coher.phaseMid);

fp2phMid = out(1).fpph.mid;
fp2ph = out(1).fpph.fp2ph;
da2ach = out(1).fpph.da2ach;
fp2phNorm = out(1).fpph.fp2phNorm;

nStates = 3; nAn = 6;

% corr.corr_achda = corr_achda; corr.lags = lags; corr.shuff_achda = shuff_achda; corr.minVal = min_val; corr.minLag = min_lag;
% coher.coher = coher_achda; coher.phase = phase_achda; coher.t = t; coher.f = f; coher.coher_shuff = coher_shuff; coher.phase_shuff = phase_shuff; coher.coherMid = coher_avg; coher.phaseMid = phase_avg;
% fpph.fp2ph = fp2ph; fpph.fp2phNorm = fp2ph_norm; fpph.da2ach = da2achph_norm; fpph.mid = mid;
% out(2).corr = corr; out(2).coher = coher; out(2).fpph = fpph;

%%
fig = figure;
fig.Position([3 4]) = [1375 1150];
clr = {'r','g','b'};

% PLOT CROSS-CORRELATION
subplot(3,3,1); hold on  
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
subplot(3,3,2); hold on
    for z = 1:nStates % iterate over behavioral states
        shadederrbar(f, nanmean(coher_achda{z},2), SEM(coher_achda{z},2), clr{z});
    end
    plot([band(1) band(1)],[0 1],'--k'); % plot lines at 0.5, 4 Hz 
    plot([band(2) band(2)],[0 1],'--k');
    shadederrbar(f, nanmean(coher_shuff{2},2), nanmean(coher_shuff{3},2) - nanmean(coher_shuff{2},2), 'k');
    legend({'imm','loc','rew','0.5Hz','4Hz'},'Location','southeast');
    xlabel('frequency');
    ylabel('coherence'); ylim([0 1]); yticks([0:0.2:1]);
    title('coherence magnitude');
    axis square; set(gca,'TickDir','out');
    
% PLOT PHASE OFFSET
subplot(3,3,3); hold on
    for z = 1:nStates % iterate over behavioral states
        shadederrbar(f, nanmean(rad2deg(phase_achda{z}),2), SEM(rad2deg(phase_achda{z}),2), clr{z}); 
    end
    plot([band(1) band(1)],[-180 180],'--k'); % plot lines at 0.5, 4 Hz
    plot([band(2) band(2)],[-180 180],'--k');
    shadederrbar(f, nanmean(rad2deg(phase_shuff{2}),2), nanmean(rad2deg(phase_shuff{3}),2) - nanmean(rad2deg(phase_shuff{2}),2), 'k');
    xlabel('frequency');
    ylabel('degrees'); ylim([-180 180]); yticks([-180:90:180]);
    title('coherence phase'); 
    axis square; set(gca,'TickDir','out');

% PLOT DA to DA phase
for y = 1:2
    subplot(3,3,3+y);
    switch y; case 1; lbl = 'ACh'; case 2; lbl = 'DA'; end
    for z = 1:nStates
        shadederrbar(fp2phMid, nanmean(fp2ph{y,z},2), SEM(fp2ph{y,z},2), clr{z});
    end
    xlabel(sprintf('%s Phase',lbl)); xlim([-180 180]); xticks([-180:90:180]);
    ylabel(sprintf('%s Fluorescence (%dF/F)',lbl)); ylim([-4 6])
    title(sprintf('%s fp norm to %s phase',lbl,lbl));
    axis square; set(gca,'TickDir','out');

end
% PLOT DA TO ACh phase
subplot(3,3,6);
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
statsb2 = {minVal, coherMid, phaseMid};
lbl = {'coefficient','coherence','phase offset'};
for x = 1:3
    subplot(3,3,6+x); hold on;
    a = statsb2{x};
    plot([0.5 3.5],[1 1],'--k');
    plot(jit,a,'.k','MarkerSize',20); % plot raw values per animal
    errorbar(nanmean(a,1), SEM(a,1),'.r', 'MarkerSize', 20);
    xlim([0.5 3.5]); xticks([1:3]); xticklabels({'imm','mov','rew'}); 
    ylabel(sprintf('%s (norm to control)',lbl{x})); ylim([0 1.5]); yticks([0:0.5:1.5]);
    p = anova1(a,[],'off');
    title(sprintf('coeff: anova p = %1.2f',p));
    axis square; set(gca,'TickDir','out');
end
    
movegui(gcf,'center');