%% LOAD
load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig4_synchrony');

%% PLOT 
x = 1; % IMMOBILITY
time = out(x).time;
nPairs = size(out(x).delta,2);

fig = figure;
fig.Position([3 4]) = [1375 800];
subplot(2,3,[1 2]); hold on
for y = 1:length(out(x).exSpikes)
    plot(out(x).exSpikes{y}, -y.*ones(length(out(x).exSpikes{y}),1), '.b', 'MarkerSize', 20);
end
ylim([-10 0]); xlabel('Time (s)');
title(sprintf('pCIN spikes %s example',out(x).lbl)); set(gca,'TickDir','out');

subplot(2,3,3); hold on
sm = 5;
shadederrbar(time, movmean(out(x).example50,sm), movmean(out(x).example95,sm), 'k');
plot(time, movmean(out(x).example,sm), 'b');
xlabel('Lag from ref. spike (s)'); xlim([-1 1]); 
ylabel('Firing rate (sp/s)'); ylim([3 7]); yticks([3:7]);
title(sprintf('CCG %s example',out(x).lbl));
axis square; set(gca,'TickDir','out');

subplot(2,3,4); hold on
[~, ii] = sort(out(x).delta(time == 0,:)); % sort in ascending order value at lag = 0
t = time(find(time == -1):find(time == 1));
h = imagesc(t, [1:nPairs], 1+out(x).delta(:,ii)', [0.6 2.2]);
colorbar; colormap(jet(256));
xlabel('Lag from ref. spike (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('pCIN pair');
h(x) = h.Parent; h(x).CLim = [0.6 2.2];
title(sprintf('CCG %s (n = %d)',out(x).lbl,nPairs));
axis square; set(gca,'TickDir','out');

subplot(2,3,5); hold on
sm = 1;
shadederrbar(time, 1+movmean(nanmean(out(x).delta50,2),sm), movmean(nanmean(out(x).delta95,2),sm), 'k');
shadederrbar(time, 1+movmean(nanmean(out(x).delta,2),sm), movmean(SEM(out(x).delta,2),sm), 'b'); 
xlabel('Lag from ref. spike (s)'); xlim([-1 1]); 
ylabel('Firing rate (norm.)'); ylim([0.9 1.4]); yticks([0:0.2:1.4]);
title(sprintf('CCG mean, max = %1.2f',max(1+nanmean(out(x).delta,2))));
axis square; set(gca,'TickDir','out');

subplot(2,3,6); hold on
ds = 2;
a = 100*sum(out(x).above95,2)/size(out(x).above95,2); 
[prop_m,prop_t] = max(a); prop_t = time(prop_t); 
a = a(1:ds:end);
b = -100*sum(out(x).below5,2)/size(out(x).below5,2); b = b(1:ds:end);
bar(time(1:ds:end), a,'FaceColor','b','FaceAlpha',0.5);
bar(time(1:ds:end), b,'FaceColor','b','FaceAlpha',0.5);
xlabel('Lag from ref. spike (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('% of CCGs'); ylim([-50 100]); xlim([-1 1]);
title(sprintf('max = %1.1f @ %d ms',prop_m, round(1000*prop_t)))
axis square; set(gca,'TickDir','out');
movegui(gcf,'center');

%% HEATMAP
figure;

[~, ii] = sort(out(x).delta(time == 0,:)); % sort in ascending order value at lag = 0
t = time(find(time == -1):find(time == 1));
nPairs = size(out(x).delta,2);
h = imagesc(t, [1:nPairs], 1+out(x).delta(:,ii)', [0.6 2.2]);
colorbar; colormap(jet(256));
xlabel('Lag from ref. spike (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('pCIN pair');
h(y) = h.Parent; h(y).CLim = [0.6 2.2];
title(sprintf('CCG %s (n = %d)',out(x).lbl,nPairs));
axis square; set(gca,'TickDir','out');