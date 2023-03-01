load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig4_fp2syncST.mat');

fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); hold on
clr = {'m','r','c','b','g'};
x = 4; cts = matExample(x).cts0; % extract counts
cts(cts > 1) = 1; 
cts = sum(cts, 1); % sum across units
time = [-1:1/50:1]';
for y = 1:length(matExample)
    iTmp = find(cts == y-1); iTmp = iTmp'; % index of coherence among units
    iTmp_imm = matExample(x).idxImm((ismember(matExample(x).idxImm, iTmp))); % index of coherence among units that occur during immobility
    pull_sta = matExample(x).sta;
    pull_sta = pull_sta - nanmean(pull_sta([1:10],:));
    shadederrbar(time, 1+nanmean(pull_sta(:,iTmp_imm),2), SEM(pull_sta(:,iTmp_imm),2), clr{y});
end
xlabel('Time to spike (s)'); xlim([-0.5 0.5]);
ylabel('ACh (norm.)');
title(sprintf('%s-#%d',matExample(x).rec,matExample(x).n),'Interpreter','none');
axis('square');

subplot(1,2,2); hold on
time = [-1:1/50:1]';
plot([0 0],[0.5 2.0],'--k'); 
shadederrbar(time, 1+nanmean(staZero,2), SEM(staZero,2), 'm'); hold on
shadederrbar(time, 1+nanmean(staMax,2), SEM(staMax,2), 'g');
legend({'1/N','N/N'});
ylabel('ACh (norm.)'); ylim([0.5 2.5]); 
xlabel('Time to spike (s)'); xlim([-0.5 0.5]);
title(sprintf('1/N vs N/N spikes (n = %d units)',size(staMax,2)));
axis('square');
movegui(fig, 'center');

%%
% above95_max = []; above95_zero = []; below5_max = []; below5_zero = [];
% for x = 1:length(mat)
%     if isnan(staZeroImm(1,x)); continue; end
%     above95_zero = [above95_zero, staZeroImm(:,x) > (staShuff50(:,x)+staShuff95(:,x))];
%     above95_max  = [above95_max,  staMaxImm(:,x)  > (staShuff50(:,x)+staShuff95(:,x))];
%     below5_zero  = [below5_zero,  staZeroImm(:,x) < (staShuff50(:,x)-staShuff95(:,x))];
%     below5_max   = [below5_max,   staMaxImm(:,x)  < (staShuff50(:,x)-staShuff95(:,x))];
% end
% ds = 2;
% above95_max = above95_max(1:ds:end,:); above95_zero = above95_zero(1:ds:end,:); 
% below5_max = below5_max(1:ds:end,:); below5_zero = below5_zero(1:ds:end,:);
% % figure; hold on
% subplot(1,2,2); hold on
% bar(time(1:ds:end), 100*sum(above95_zero,2)/size(above95_zero,2),'FaceColor','m','FaceAlpha',0.5);
% bar(time(1:ds:end), -100*sum(below5_zero,2)/size(below5_zero,2),'FaceColor','m','FaceAlpha',0.5);
% bar(time(1:ds:end), 100*sum(above95_max,2)/size(above95_max,2),'FaceColor','g','FaceAlpha',0.5);
% bar(time(1:ds:end), -100*sum(below5_max,2)/size(below5_max,2),'FaceColor','g','FaceAlpha',0.5);
% xlabel('Time to spike (s)'); xlim([-1 1]);
% ylabel('% of units'); ylim([-100 100]);
% title(sprintf('% of units (n = %d)',size(below5_max,2)))
% axis('square');