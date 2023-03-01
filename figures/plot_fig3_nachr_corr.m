load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig3_nachr_beh','cannula');

%% PLOT STATISTICS
fig = figure;
fig.Position(4) = 1150;
clr = {'r','g','b'};
jit = []; % jitter x-values for plotting
for z = 1:nStates
    j1 = z-0.25; j2 = z+0.25; jit(:,z) = j1 + (j2-j1).*rand(nAn,1); % jitter x-values for plotting
end
minVal = cannula(2).corrMinVal./cannula(1).corrMinVal; % Normalize to saline
coherMid = cannula(2).coherMid./cannula(1).coherMid; % Normalize to saline
phaseMid = cannula(2).phaseMid./cannula(1).phaseMid; % Normalize to saline
statsCan = {minVal, coherMid, phaseMid};
lbl = {'coefficient','coherence','phase offset'};
for x = 1:3
    subplot(3,1,x); hold on;
    a = statsCan{x};
    plot([0.5 3.5],[1 1],'--k');
    plot(jit,a,'.k','MarkerSize',20); % plot raw values per animal
    errorbar(nanmean(a,1), SEM(a,1),'.r', 'MarkerSize', 20);
    xlim([0.5 3.5]); xticks([1:3]); xticklabels({'imm','mov','rew'}); 
    ylabel(sprintf('%s (norm to saline)',lbl{x})); ylim([0 1.5]); yticks([0:0.5:1.5]);
    p = anova1(a,[],'off');
    title(sprintf('coeff: anova p = %1.2f',p));
    axis square; set(gca,'TickDir','out');
end
movegui(gcf,'center');