load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig3_b2floxTraces.mat');

a = [2 2 1];
b = [1620 1630; 1580 1590; 815 825];
fig = figure;
fig.Position(3) = 1375;
for y = 1:3
    x = a(y);
    fp_mat = [beh(x).FP{1}, beh(x).FP{2}]; % extract photometry signal from structure
    fp_mat = fp_mat - nanmean(fp_mat); % subtract mean of entire photometry signal, now centered
    acc = getAcc(beh(x).vel); % acceleration
    sp(y) = subplot(1,3,y); hold on
        plot(beh(x).time, fp_mat(:,1), 'g'); 
        plot(beh(x).time, fp_mat(:,2), 'm');
        xlabel('Time (s)'); xlim(b(y,:));
        ylabel('Photometry (%dF/F)'); ylim([-6 10])
        plot(beh(x).time, acc - 5, 'k');
        title('b2flox');
        axis square; set(gca,'TickDir','out');
end
linkaxes(sp,'y');
movegui(gcf,'center');
