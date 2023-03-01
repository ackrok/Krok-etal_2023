load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig4_tuning');

fig = figure; fig.Position(3) = 1000;
sm = 10;
subplot(1,2,1); hold on 
phasedistros = movmean([mat.phasedistros],sm); % Spike probability to phase
phasedistros(:,[mat.p] >= 0.05) = []; % Remove untuned units
phasedistros = normalize(phasedistros,'range'); % Normalize
plot(rad2deg(phasebins), cinExample, 'c');
shadederrbar(rad2deg(phasebins), nanmean(phasedistros,2), SEM(phasedistros,2), 'b');
ylabel('Spike probability (norm.)'); ylim([0 1]); yticks([0:0.5:1]);
ylim([0 1]);
yyaxis right
plot(rad2deg(phasebins),  cos(phasebins), '--k');    
title(sprintf('n = %d units / %d are tuned',size(phasedistros,2),size([mat.phasedistros],2)));
xlabel('LFP phase'); xlim([0 360]); xticks([0:90:360]);
axis('square');

subplot(1,2,2); % Phase preferences for all units
pref = [mat.pref]; pref([mat.p] >= 0.05) = []; % Remove untuned units 
rose(deg2rad(pref));
movegui(gcf,'center');