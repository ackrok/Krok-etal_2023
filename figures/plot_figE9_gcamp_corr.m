load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_gcamp-DLS-SNc');

%% inputs
winPlot = [-1 1]; % window for plotting

%% COMPUTE CROSS_CORRELATION
mat = struct;
for x = 1:length(beh); y = [2 1]; %CHANGE - which FP signal to use as reference
    %% extract signals
    Fs = beh(x).Fs; 
    fp_mat = [];
    fp_mat(:,1) = beh(x).FP{y(1)}; % extract photometry signal from structure
    fp_mat(:,2) = beh(x).FP{y(2)};
    fp_mat = fp_mat - nanmean(fp_mat); % center traces
    idx_imm = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    fp_sub = fp_mat(idx_imm,:); % extract signal during immobility
    [c,lags] = xcorr(fp_sub(:,1), fp_sub(:,2), 10*Fs, 'coeff'); % cross-correlation during immobility
    tmp_shuff = []; 
    for s = 1:50
        fp_shift = fp_sub(randperm(length(fp_sub)),2); % shuffle photometry values
        tmp_shuff(:,s) = xcorr(fp_sub(:,1), fp_shift, 10*Fs, 'coeff'); % cross-correlation of shuffled DLS signal
    end
    mat(x).rec = beh(x).rec;
    mat(x).ref = beh(x).FPnames{y(1)};
    mat(x).run = beh(x).FPnames{y(2)};
    mat(x).lags = lags./Fs;
    mat(x).c = c;
    mat(x).shuff = prctile(tmp_shuff, [2.5 50 97.5], 2);
end
fprintf('Cross-correlation GCaMP6f DLS/SNc DONE!\n');

%% EXTRACT MAXIMUM COEFFICIENT by mouse
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); nAn = length(uni);
mat_corr = []; % initialize output
mat_shuff5 = []; mat_shuff50 = [];
for x = 1:nAn
    idx = find(strcmp(rec,uni{x})); % identify all recordings from this animal
    [~,idx2] = max(max([mat(idx).c]));
    mat_corr(:,x) = [mat(idx(idx2)).c]; % extract recording
    mat_corr(:,x) = mat_corr(:,x) - nanmean(mat_corr(find([lags/Fs] == -2):find([lags/Fs] == -1),x));
    mat_shuff5(:,x) = mat(idx(idx2)).shuff(:,1);
    mat_shuff50(:,x) = mat(idx(idx2)).shuff(:,2);
end
[mm, ii] = max(mat_corr); % peak correlation coefficient
ii = lags(ii)/Fs*1000; % convert to milliseconds

%% PLOT
fig = figure; fig.Position(3) = 1375;
subplot(1,3,1); hold on
x = 9; % example recording
plot(lags/Fs, mat(x).c, 'm');
shadederrbar(lags/Fs, mat(x).shuff(:,2), mat(x).shuff(:,1), 'k');
plot([0 0],[-0.2 1],'--k');
title('example corr imm DLS/SNc'); axis square
xlabel('Latency to SNc (s)'); xlim(winPlot); xticks([winPlot(1):0.5:winPlot(2)]);
ylabel('Coefficient'); ylim([-0.2 0.6]); yticks([-0.2:0.2:1]);

subplot(1,3,2); hold on
shadederrbar(lags/Fs, nanmean(mat_corr,2), SEM(mat_corr,2), 'm');
shadederrbar(lags/Fs, nanmean(mat_shuff50,2), nanmean(mat_shuff5,2), 'k');
plot([0 0],[-0.2 1],'--k');
title('average corr imm DLS/SNc'); axis square
xlabel('Latency to SNc (s)'); xlim(winPlot); xticks([winPlot(1):0.5:winPlot(2)]);
ylabel('Coefficient'); ylim([-0.2 0.6]); yticks([-0.2:0.2:1]);

subplot(1,3,3); hold on
clr = hsv(nAn);
for x = 1:nAn
    plot(ii(x)/1000, mm(x), '.', 'MarkerSize', 20, 'Color', clr(x,:)); % plot individual data points
end
plot([0 0],[0 1],'--k'); % vertical line at lag = 0
errorbar(nanmean(ii)/1000, nanmean(mm), ...
    SEM(mm,2), SEM(mm,2), ...
    SEM(ii,2)/1000, SEM(ii,2)/1000,...
    '.k', 'MarkerSize', 20); % error bar with SEM
title(sprintf('Peak: %1.2f +/- %1.2f || Lag: %1.2f +/- %1.2f ms', nanmean(mm), nanstd(mm), nanmean(ii), nanstd(ii))); 
xlabel('Latency to SNc (s)'); xlim(winPlot); xticks([winPlot(1):0.5:winPlot(2)]);
ylabel('Coefficient'); ylim([0 1]); yticks([0:0.2:1]);
legend(uni);
axis square
movegui(gcf,'center');

%% TESTING EXAMPLE
figure;
for x = 1:length(mat)
    subplot(3,4,x); hold on
    plot([0 0],[-0.2 1],'--k'); plot([-1 1],[0 0],'--k');
    plot(lags/Fs, mat(x).c); xlim([-1 1]); ylim([-0.2 1]);
    title(sprintf('%d: r = %1.2f',x,max(mat(x).c)));
end

%% PLOT TRACE
x = 3;

figure;
sp(1) = subplot(3,1,1); 
plot(beh(x).time, beh(x).FP{1}, 'b'); 
title(sprintf('%s - %s',beh(x).rec,beh(x).FPnames{1}));
% ylim([-5 20]);
sp(2) = subplot(3,1,2);
plot(beh(x).time, beh(x).FP{2}, 'm'); 
title(sprintf('%s - %s',beh(x).rec,beh(x).FPnames{2}));
% ylim([-5 60]);
sp(3) = subplot(3,1,3);
plot(beh(x).time, getAcc(beh(x).vel), 'k');
ylim([-1 1]);
linkaxes(sp,'x');
% xlim([4.5 14.5]);
