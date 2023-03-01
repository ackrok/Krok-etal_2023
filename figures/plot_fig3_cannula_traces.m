load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig3_cannulaTraces.mat');

% saline segments: [2185 2195], [1922.5 1932.5], [1560 1570])
% d1d2R segments: [1191.5 1199.5]
% iGluR segments: [3125 3135]
% nAChR segments: [2555 2565], [3405 3415]

%% Saline, D1/2R antagonist
fig = figure; fig.Position(3) = 1000;
a = [1 2];
b = [1922.5 1932.5; 1191.5 1199.5];
clr = {'k','b'};
for c = 1:2; x = a(c);
    subplot(2,length(a),c); hold on
    plot(beh(x).time, beh(x).FP{1}, clr{c}); % plot ACh
    title(sprintf('%s',beh(x).inf)); set(gca,'TickDir','out');
    ylabel('(%dF/F)'); ylim([-5 22])
    xlim(b(c,:));
    
    subplot(2,length(a),length(a)+c);
    plot(beh(x).time, getAcc(beh(x).vel), 'Color', [0.5 0.5 0.5]); % plot acceleration
    title(sprintf('%s',beh(x).inf)); set(gca,'TickDir','out');
    ylabel('Acceleration'); ylim([-1 1]);
    xlim(b(c,:));
end
movegui(gcf,'center');

%% Saline, iGluR antagonist
fig = figure; fig.Position(3) = 1000;
a = [1 3];
b = [2185 2195; 3125 3135];
for c = 1:2; x = a(c);
    subplot(2,length(a),c); hold on
    plot(beh(x).time, beh(x).FP{1}, clr{c}); % plot ACh
    title(sprintf('%s',beh(x).inf)); set(gca,'TickDir','out');
    ylabel('(%dF/F)'); ylim([-5 22])
    xlim(b(c,:));
    
    subplot(2,length(a),length(a)+c);
    plot(beh(x).time, getAcc(beh(x).vel), 'Color', [0.5 0.5 0.5]); % plot acceleration
    title(sprintf('%s',beh(x).inf)); set(gca,'TickDir','out');
    ylabel('Acceleration'); ylim([-1 1]);
    xlim(b(c,:));
end
movegui(gcf,'center');