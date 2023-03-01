load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_figSX_gcamp_fft_data.mat')
%%
f = out(1).f; flog = out(1).flog; range_auc = out(1).range_auc;
sub_comp = cell(1,2); sub_comp{1} = out(1).sub_comp; sub_comp{2} = out(2).sub_comp; 
auc = []; auc(:,1) = out(1).auc(:); auc(:,2) = out(2).auc(:);
lbl = {out.lbl_fp};
ds = 50; % Downsample when plotting so figure is smaller in size
nAn = size(sub_comp{1},2); % Number of recordings and/or animals

fig = figure; fig.Position(3) = 1375;
for x = 1:2
    subplot(1,3,x); hold on
    plot(flog(1:ds:end), sub_comp{x}((1:ds:end),:)); % Plot FFT for each recording
    plot([-2 2],[0 0],'-','Color',[0 0 0 0.3]); % Line at 0 power
    plot([flog(f == range_auc(1)) flog(f == range_auc(1))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 0.5 Hz
    plot([flog(f == range_auc(2)) flog(f == range_auc(2))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 4 Hz 
    xlabel('Frequency'); ylabel('Power (a.u.)');
    xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
    title(sprintf('%s FFT immobility',lbl{x})); axis square
end
subplot(1,3,3); hold on
    clr = {'k','b'};
    for x = 1:2
        j1 = 0.8; j2 = 1.2; jit = j1 + (j2-j1).*rand(nAn,1); % jitter
        plot(jit, auc(:,x), '.', 'Color', clr{x}, 'MarkerSize', 20); % Plot area under curve power data points for each recording
        errorbar(nanmean(auc(:,x)), SEM(auc(:,x),1), '.', 'Color', clr{x}, 'MarkerSize', 20); % Plot area under curve power averaged
    end
    xlim([0.5 1.5]); xticks([1]); xticklabels({'AUC'});
    ylabel('Power (a.u.)'); ylim([0 0.5]);
    title(sprintf('AUC %1.1f-%dHz - %s: (%1.2f) vs %s (%1.2f)',range_auc(1),range_auc(2),lbl{1},nanmean(auc(:,1)),lbl{2},nanmean(auc(:,2)))); axis square
movegui(gcf,'center');