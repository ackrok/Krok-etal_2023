%%
load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_figS3_dms');

%% ACh, DA power comparing DLS and DMS antagonist
f = outFft(1).f; flog = outFft(1).flog; % extract frequency domain vectors
ds = 50; % downsample vectors for plotting
range_auc = [0.5 4]; % Range for calculation of area under the curve, in Hz
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_auc = [find(f_sub == range_auc(1)):find(f_sub == range_auc(2))]; % extract AUC from [0.5 4] Hz
clr = {'k','m';'k','g'};
fig = figure; fig.Position([3 4]) = [850 750];
for x = 1:length(outFft)
    plotme = {outFft(x).dls, outFft(x).dms};
    sub_mat = cell(1,length(plotme));
    auc = [];
    for z = 1:length(plotme)
        for y = 1:size(plotme{z},2)
            sub_mat{z}(:,y) = plotme{z}(:,y) - nanmean(outFft(x).fluor,2); % subtract FFT of stable fluorophore
            auc(y,z) = trapz(sub_mat{z}(r_auc,y))/length(r_auc);
        end
    end
    auc(auc == 0) = nan; auc(auc < 0) = 0;
    %
    subplot(2,2,x); hold on
    for z = 1:length(sub_mat)
        shadederrbar(flog(1:ds:end), nanmean(sub_mat{z}(1:ds:end,:),2), SEM(sub_mat{z}(1:ds:end,:),2), clr{x,z}); % plot average FFT outAntput
    end
        plot([-2 2],[0 0],'-','Color',[0 0 0 0.3]); % Line at 0 power
        plot([flog(f == range_auc(1)) flog(f == range_auc(1))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 0.5 Hz
        plot([flog(f == range_auc(2)) flog(f == range_auc(2))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 4 Hz 
        xlabel('Frequency'); ylabel('Power (a.u.)');
        xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
        legend(outFft(x).lbl);
        title(sprintf('%s: %s vs %s',outFft(x).fplbl,outFft(x).lbl{1},outFft(x).lbl{2}));
        axis square
    %
    [~,p] = ttest(auc(:,1),auc(:,2));
    subplot(2,2,x+2); hold on
        plot([1;2].*ones(2,size(auc,1)), auc', '.-', 'MarkerSize', 20, 'Color', clr{x,z}); % Plot area under curve power data points for each recording
        errorbar([0.75 2.25], nanmean(auc), SEM(auc,1), '.k', 'MarkerSize', 20);
        xlim([0.5 2.5]); xticks([1:3]); xticklabels(outFft(x).lbl);
        ylabel('Power (a.u.)'); ylim([0 0.5]); yticks([0:0.1:0.5]);
        title(sprintf('p = %1.4f (n = %d)',p,size(auc,1)));
        axis square
end
movegui(gcf,'center');
