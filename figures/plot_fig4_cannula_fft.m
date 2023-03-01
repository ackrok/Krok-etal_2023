% load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig4_iglur2_fft_data');
% load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig3_nachr_fft_data');
% load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_figX_machr_fft_data');
% 
% % out = struct;
% % out(1).lbl_fp = 'ACh'; out(1).sub_comp = sub1ach; out(1).lbl_inf = lbl_inf;
% % out(2).lbl_fp = 'DA'; out(2).sub_comp = sub1da; out(2).lbl_inf = lbl_inf;
% % out(1).f = f; out(1).flog = flog; 
% % out(1).range_auc = [0.5 4]; out(2).range_auc = [0.5 4];

%% PLOT
for z = 1:2
    lbl_fp = out(z).lbl_fp;
    sub_comp = out(z).sub_comp;
    lbl_inf = out(z).lbl_inf;
    auc = out(z).auc; nAn = size(auc,1);
    f = out(1).f; flog = out(1).flog; range_auc = out(1).range_auc;
    
    ds = 50; % Downsample when plotting so figure is smaller in size
    fig = figure; fig.Position(3) = 1000;
    switch lbl_fp
        case 'ACh'; clr = {'k','g','r','m','b','c'}; 
            clr2 = [0 0 0 0.2; 0.05 0.75 0.45 0.2; 0 0 1 0.2; 0.05 0.75 0.45 0.2; 0 0 1 0.2];
        case 'DA'; clr = {'k','m','r','g','b','c'}; 
            clr2 = [0 0 0 0.2; 1 0 1 0.2; 0 0 1 0.2; 0.05 0.75 0.45 0.2; 0 0 1 0.2];
    end
    subplot(1,2,1); hold on
        for y = 1:length(sub_comp)
            plotme = sub_comp{y}((1:ds:end),:);
            shadederrbar(flog(1:ds:end), nanmean(plotme,2), SEM(plotme,2), clr{y}); % Plot FFT averaged across all recordings
        end
        plot([-2 2],[0 0],'-','Color',[0 0 0 0.3]); % Line at 0 power
        plot([flog(f == range_auc(1)) flog(f == range_auc(1))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 0.5 Hz
        plot([flog(f == range_auc(2)) flog(f == range_auc(2))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 4 Hz 
        ylabel('Power (a.u.)'); yticks([-0.1:0.1:0.4]);
        xlabel('Frequency'); xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
        title(sprintf('FFT (n = %d): %s', nAn, lbl_fp)); axis square
        legend(lbl_inf);
    subplot(1,2,2); hold on
        plot(auc', '--.k', 'MarkerSize', 20); 
        errorbar(-0.25+1:size(auc,2),nanmean(auc), SEM(auc,1), '.', 'MarkerSize', 20, 'Color', clr{2});
        xlim([0.5 0.5+size(auc,2)]); xticks([1:size(auc,2)]); xticklabels(lbl_inf);
        ylabel('Power (a.u.)'); ylim([0 0.5]); yticks([0:0.1:0.5]);
        [~,p] = ttest(auc(:,1),auc(:,2));
        title(sprintf('AUC p = %1.3f',p)); axis square
    movegui(gcf,'center');
end