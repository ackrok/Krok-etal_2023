males = [1 2 3 6 8 9]; females = [4 5 7 10];

load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig1_fft')

%%
f = out(1).f; flog = out(1).flog; % extract frequency domain vectors
ds = 50; % downsample vectors for plotting
range_auc = [0.5 4]; % Range for calculation of area under the curve, in Hz
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_auc = [find(f_sub == range_auc(1)):find(f_sub == range_auc(2))]; % extract AUC from [0.5 4] Hz
clr = {'b','r'};
fig = figure; fig.Position(3) = 1375;
for x = 1:2
    sub_mat = [];
    for y = 1:size(out(x).control,2)
        sub_mat(:,y) = out(x).control(:,y) - nanmean(out(x).fluor,2); % subtract FFT of stable fluorophore
        auc{x}(y) = trapz(sub_mat(r_auc,y))/length(r_auc);
    end
    % [~,ii] = max(nanmean(sub_mat(:,males),2)); mode1 = f(ii);
    % [~,ii] = max(nanmean(sub_mat(:,females),2)); mode2 = f(ii);
    [~,ii] = max(sub_mat(:,males)); mode1(1) = f(round(nanmean(ii))); mode1(2) = f(round(SEM(ii,2)));
    [~,ii] = max(sub_mat(:,females)); mode2(1) = f(round(nanmean(ii))); mode2(2) = f(round(SEM(ii,2)));
    subplot(1,3,x); hold on
        shadederrbar(flog(1:ds:end), nanmean(sub_mat(1:ds:end,males),2), SEM(sub_mat(1:ds:end,males),2), clr{1});
        shadederrbar(flog(1:ds:end), nanmean(sub_mat(1:ds:end,females),2), SEM(sub_mat(1:ds:end,females),2), clr{2}); % plot average FFT output
        plot([-2 2],[0 0],'-','Color',[0 0 0 0.3]); % Line at 0 power
        plot([flog(f == range_auc(1)) flog(f == range_auc(1))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 0.5 Hz
        plot([flog(f == range_auc(2)) flog(f == range_auc(2))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 4 Hz 
        xlabel('Frequency'); ylabel('Power (a.u.)');
        xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
        axis square
        title(sprintf('%s (M: %1.1f +/- %1.1f Hz, F: %1.1f +/- %1.1f Hz)',out(x).fplbl,mode1(1),mode1(2),mode2(1),mode2(2)));
        legend({'males','females'});
        
end

subplot(1,3,3); hold on
clr = {'b','r','b','r'}; 
auc2 = cell(1,4); 
auc2{1} = auc{1}(males); auc2{2} = auc{1}(females); 
auc2{3} = auc{2}(males); auc2{4} = auc{2}(females);
for z = 1:length(auc2)
    j1 = -0.2+z; j2 = 0.2+z; jit = j1 + (j2-j1).*rand(length(auc2{z}),1); % jitter
    plot(jit, auc2{z}', '.', 'MarkerSize', 20, 'Color', clr{z}); % Plot area under curve power data points for each recording
end
    errorbar([1:length(auc2)], cellfun(@nanmean, auc2), cellfun(@nanstd, auc2)./sqrt(cellfun(@length, auc2)), '.k', 'MarkerSize', 20);
    xlim([0.5 4.5]); xticks([1:4]); xticklabels({'ACh M','ACh F','DA M','DA F'});
    ylabel('Power (a.u.)');
    axis square
    [~,p] = ttest2(auc2{1},auc2{2});[~,p(2)] = ttest2(auc2{3},auc2{4});
    title(sprintf('ACh M/F p = %1.2f, DA M/F p = %1.2f',p(1),p(2)));
