load('C:\Users\Anya\Desktop\FP_LOCAL\krok_fig1_fft')

%% Figure 1I, 1M
f = out(1).f; flog = out(1).flog; % extract frequency domain vectors
ds = 50; % downsample vectors for plotting
clr = {'g','m'};
for x = 1:length(out)
    plotme = {out(x).fluor, out(x).control, out(x).antag, out(x).lesion};
    lbl = {'fluorophore','control','antagonist','lesion'};
    fig = figure; fig.Position(3) = 450;
    for y = 1:4
    subplot(2,2,y);
        shadederrbar(flog(1:ds:end), nanmean(plotme{y}(1:ds:end,:),2), SEM(plotme{y}(1:ds:end,:),2), clr{x}); % plot average FFT output for each condition
        xlabel('Frequency'); ylabel('Power (a.u.)');
        xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
        title(sprintf('%s',lbl{y})); axis square
    end
    subplot(2,2,2); hold on;
    plot(flog(1:ds:end), nanmean(plotme{1}(1:ds:end,:),2), '--k'); % add stable fluorophore FFT output to control subplot
end

%% Figure 1J, K, N, O
f = out(1).f; flog = out(1).flog; % extract frequency domain vectors
ds = 50; % downsample vectors for plotting
range_auc = [0.5 4]; % Range for calculation of area under the curve, in Hz
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_auc = [find(f_sub == range_auc(1)):find(f_sub == range_auc(2))]; % extract AUC from [0.5 4] Hz
clr = {'g','m'};
figure;
for x = 1:length(out)
    sub_mat = [];
    for y = 1:size(out(x).control,2)
        sub_mat(:,y) = out(x).control(:,y) - nanmean(out(x).fluor,2); % subtract FFT of stable fluorophore
    end
    subplot(2,2,x); hold on
        shadederrbar(flog(1:ds:end), nanmean(sub_mat(1:ds:end,:),2), SEM(sub_mat(1:ds:end,:),2), clr{x}); % plot average FFT output
        plot([-2 2],[0 0],'-','Color',[0 0 0 0.3]); % Line at 0 power
        plot([flog(f == range_auc(1)) flog(f == range_auc(1))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 0.5 Hz
        plot([flog(f == range_auc(2)) flog(f == range_auc(2))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 4 Hz 
        xlabel('Frequency'); ylabel('Power (a.u.)');
        xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
        axis square
    plotme = {out(x).control, out(x).antag, out(x).lesion};
    auc = cell(length(plotme),1);
    subplot(2,2,x+2); hold on
    for z = 1:length(plotme)
        for y = 1:size(plotme{z},2)
            sub_mat = plotme{z}(:,y) - nanmean(out(x).fluor,2); % subtract FFT of stable fluorophore
            auc{z}(y) = trapz(sub_mat(r_auc))/length(r_auc);
        end
        auc{z}(auc{z} == 0) = nan; 
        auc{z}(auc{z} < 0) = 0;
        j1 = -0.2+z; j2 = 0.2+z; jit = j1 + (j2-j1).*rand(length(auc{z}),1); % jitter
        plot(jit, auc{z}, '.', 'MarkerSize', 20, 'Color', clr{x}); % Plot area under curve power data points for each recording
    end
    errorbar([1:length(plotme)], cellfun(@nanmean, auc), cellfun(@nanstd, auc)./sqrt(cellfun(@length, auc)), '.k', 'MarkerSize', 20);
    xlim([0.5 3.5]); xticks([1:3]); xticklabels({'control','antag','lesion'});
    ylabel('Power (a.u.)');
    axis square
end
