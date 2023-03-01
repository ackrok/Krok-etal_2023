load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_figOF_corr2.mat');
load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_figOF_coher2.mat');
load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_figOF_fp2ph2.mat');

% load('C:\Users\Anya\Desktop\FP_LOCAL\openField\beh_OF_SummaryNov22.mat');
% plot_achda_2
% uni = {'AK231','AK232','AK233','AK237'};
% save('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_figOF_corr2.mat','lags','corr_achda','shuff_achda','uni');
% save('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_figOF_coher2.mat','f','t','coher_achda','coher_shuff','phase_achda','phase_shuff','uni');
% save('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_figOF_fp2ph2.mat','mid','fp2ph','fp2ph_norm','da2achph_norm','uni');

%% CROSS_CORRELATION
fig = figure; fig.Position(3) = 1375;
clr = {'r','g','b'};
Fs = 50;
nStates = 2;
nAn = length(uni);
% PLOT CROSS-CORRELATION
subplot(1,3,1); hold on; 
    plot([0 0],[-1 0.5],'--k');
    shadederrbar(lags/Fs, nanmean(shuff_achda{1,2},2) - nanmean(nanmean(shuff_achda{1,2},2)),...
        nanmean(shuff_achda{1,1},2), 'k'); hold on
    for z = 1:nStates
        shadederrbar(lags/Fs, nanmean(corr_achda{z,1},2), SEM(corr_achda{z,1},2), clr{z});
    end
    xlabel('lag from DA (s)'); xlim([-1 1]);
    ylabel('coefficient'); ylim([-1 0.5]); yticks([-1:0.25:1]);
    title(sprintf('DA/ACh corr (n = %d mice)',nAn));
    axis square
% STATISTICS: ACh3.0 and rDA photometry are mostly anti-correlated with a negative lag
min_val = []; min_lag = []; % initialize output
for z = 1:nStates % iterate over behavioral states
    [min_val(:,z), ii] = min(corr_achda{z,1}); % cross-correlation maximal value for each animal during each behavioral state
    min_lag(:,z) = lags(ii)./Fs; % latency at minimum, in seconds
end
jit = []; % jitter x-values for plotting
for z = 1:nStates
    j1 = z-0.25; j2 = z+0.25; jit(:,z) = j1 + (j2-j1).*rand(nAn,1); % jitter x-values for plotting
end
% PLOT STATISTICS - COEFFICIENT
subplot(1,3,2); hold on;
    plotme = min_val;
    plot(jit,plotme,'.k','MarkerSize',20); % plot raw values per animal
    errorbar(nanmean(plotme,1), SEM(plotme,1),'.r', 'MarkerSize', 20);
    xlim([0.5 2.5]); xticks([1:2]); xticklabels({'imm','mov'}); 
    ylabel('coefficient'); ylim([-1 0]); yticks([-1:0.2:1]);
    [~,p] = ttest(plotme(:,1),plotme(:,2));
    title(sprintf('coeff: ttest p = %1.3f',p));
    axis square
% PLOT STATISTICS - LATENCY    
subplot(1,3,3); hold on;
    min_lag(abs(min_lag) == 5) = nan; % if lag is +/-5s, then adjust to be nan
    plotme = min_lag.*1000; % adjust lag at minimum to be in milliseconds, from seconds
    plot(jit,plotme,'.k','MarkerSize',20); % plot raw values per animal
    errorbar(nanmean(plotme,1), SEM(plotme,1),'.r', 'MarkerSize', 20);
    xlim([0.5 2.5]); xticks([1:2]); xticklabels({'imm','mov'}); 
    ylabel('lag (ms)'); ylim([-250 0]); yticks([-250:50:0]);
    [~,p] = ttest(plotme(:,1),plotme(:,2));
    title(sprintf('lag: ttest p = %1.3f',p));
    axis square
movegui(gcf,'center');

%% COHERENCE PHASE
fig = figure; fig.Position([3 4]) = [1000 860];
clr = {'r','g','b'};
band = [0.5 4]; % frequency band: 0.5-4 Hz
% PLOT COHERENCE
% plot coherence magnitude, averaged across all animals for each behavioral state
subplot(2,2,1); hold on
    for z = 1:nStates % iterate over behavioral states
        shadederrbar(f, nanmean(coher_achda{z},2), SEM(coher_achda{z},2), clr{z});
    end
    plot([band(1) band(1)],[0 1],'--k'); % plot lines at 0.5, 4 Hz 
    plot([band(2) band(2)],[0 1],'--k');
    shadederrbar(f, nanmean(coher_shuff{2},2), nanmean(coher_shuff{3},2) - nanmean(coher_shuff{2},2), 'k');
    legend({'imm','loc','0.5Hz','4Hz'});
    xlabel('frequency');
    ylabel('coherence'); ylim([0 1]); yticks([0:0.2:1]);
    title(sprintf('DA/ACh coherence (n = %d mice)',nAn));
    axis square
% PLOT PHASE OFFSET
% plot phase (in degrees), averaged across all animals for each behavioral state
subplot(2,2,2); hold on
    for z = 1:nStates % iterate over behavioral states
        shadederrbar(f, nanmean(rad2deg(phase_achda{z}),2), SEM(rad2deg(phase_achda{z}),2), clr{z}); 
    end
    plot([band(1) band(1)],[-180 180],'--k'); % plot lines at 0.5, 4 Hz
    plot([band(2) band(2)],[-180 180],'--k');
    shadederrbar(f, nanmean(rad2deg(phase_shuff{2}),2), nanmean(rad2deg(phase_shuff{3}),2) - nanmean(rad2deg(phase_shuff{2}),2), 'k');
    % legend({'imm','','loc','','rew','','0.5Hz','4Hz'});
    xlabel('frequency');
    ylabel('degrees'); ylim([-180 180]); yticks([-180:90:180]);
    title('coherence phase'); 
    axis square
% STATISTICS
r = [6:42]; % range for 0.5-4Hz
coher_avg = []; phase_avg = []; % initialize output
for z = 1:nStates % iterate over behavioral states
    coher_avg(:,z) = median(coher_achda{z}(r,:)); % median within frequency band
    phase_avg(:,z) = median(phase_achda{z}(r,:)); % median within frequency band
end
jit = []; % jitter x-values for plotting
for z = 1:nStates
    j1 = z-0.25; j2 = z+0.25; jit(:,z) = j1 + (j2-j1).*rand(nAn,1); % jitter x-values for plotting
end
% PLOT STATISTICS - COHERENCE MAGNITUDE
subplot(2,2,3); hold on
    plotme = coher_avg;
    plot(jit,plotme,'.k','MarkerSize',20); % plot raw values per animal
    errorbar(nanmean(plotme,1), SEM(plotme,1),'.r', 'MarkerSize', 20);
    xlim([0.5 2.5]); xticks([1:3]); xticklabels({'imm','mov','rew'}); 
    ylabel('coherence'); ylim([0 1]); yticks([0:0.2:1]);
    [~,p] = ttest(plotme(:,1),plotme(:,2));
    title(sprintf('coherence: ttest p = %1.2f',p));
    axis square
% PLOT STATISTICS - PHASE OFFSET
subplot(2,2,4); hold on
    plotme = rad2deg(phase_avg);
    plot(jit,plotme,'.k','MarkerSize',20); % plot raw values per animal
    errorbar(nanmean(plotme,1), SEM(plotme,1),'.r', 'MarkerSize', 20);
    xlim([0.5 2.5]); xticks([1:3]); xticklabels({'imm','mov','rew'}); 
    ylabel('phase (deg)'); ylim([-180 0]); yticks([-180:45:180]);
    [~,p] = ttest(plotme(:,1),plotme(:,2));
    title(sprintf('phase: ttest p = %1.2f',p));
    axis square
%
movegui(gcf,'center');
%
%% PHOTOMETRY TO PHASE
fig = figure; 
fig.Position([3 4]) = [1375 800];
clr = {'r','g','b'};
% PLOTTING
for y = 1:2
    switch y; case 1; lbl = 'ACh'; case 2; lbl = 'DA'; end
    plotme = y; sp(y) = subplot(2,3,plotme); hold on
    for z = 1:nStates
        shadederrbar( mid, nanmean(fp2ph{y,z},2), SEM(fp2ph{y,z},2), clr{z});
    end
    xlabel(sprintf('%s Phase',lbl)); xlim([-180 180]); xticks([-180:90:180]);
    ylabel(sprintf('%s Fluorescence (dF/F)',lbl));
    legend({'imm','mov','rew'});
    title(sprintf('%s fp to %s phase (n = %d)',lbl,lbl,size(fp2ph{1,1},2))); axis square

    plotme = y+2; subplot(2,3,plotme); hold on
    for z = 1:nStates
        shadederrbar( mid, nanmean(fp2ph_norm{y,z},2), SEM(fp2ph_norm{y,z},2), clr{z});
    end
    xlabel(sprintf('%s Phase',lbl)); xlim([-180 180]); xticks([-180:90:180]);
    ylabel(sprintf('%s Fluorescence (norm)',lbl)); ylim([0 1]);
    title(sprintf('%s fp norm to %s phase',lbl,lbl)); axis square

    plotme = y+4; subplot(2,3,plotme); hold on
    sm = 5;
    for z = 1:nStates
        shadederrbar( mid, movmean(nanmean(da2achph_norm{y,z},2),sm), movmean(SEM(da2achph_norm{y,z},2),sm), clr{z});
    end
    switch y; case 1; xx = 2; case 2; xx = 1; end
    plot( mid, nanmean(fp2ph_norm{xx,1},2), '--k');
    xlabel(sprintf('%s Phase',lbl)); xlim([-180 180]); xticks([-180:90:180]);
    ylabel('Fluorescence (norm)'); ylim([0 1]);
    title(sprintf('fp norm to %s phase',lbl)); axis square
end
movegui(gcf, 'center')