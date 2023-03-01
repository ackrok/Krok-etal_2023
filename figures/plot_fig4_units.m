load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig4_units')

%% Scatterplot of spike properties, all units
spkDuration = [units.spkDur]; % extract spike duration for each unit
firingRate    = [units.fr];  % extract firing rate for each unit
phasicIndex = [units.phasicIndex]; % extract phasic activity index for each unit
lbl = {'pSPN','pCIN','other'}; % unit groups
clr = [0.5 0.5 0.5; 0 0 1; 0 1 0]; 
figure; hold on
for x = unique([units.label])
    idx = find([units.label] == x); % find index for all units of one group
    plot3(phasicIndex(idx), spkDuration(idx), firingRate(idx),... % plot 3D data point for each unit
        '.','MarkerSize',10,'Color',clr(x,:),'DisplayName', lbl{x});
end
xlabel('Phasic activity index'); xticks([0:0.25:1]); xlim([0 1]);
ylabel('Spike duration (ms)'); yticks([0:1:3]); ylim([0 3]); 
zlabel('Firing rate (sp/s)'); zticks([0:10:50]); zlim([0 50]);
legend; grid on; view(45,45);

%% Average waveforms
figure; hold on
t = [-4:1/30:4]; % time vector, in ms, for plotting
lbl = {'pSPN','pCIN','other'}; % unit groups
clr = {'k','b','g'}; 
for x = unique([units.label])
    sub = units(find([units.label] == x)); % sub-structure
    wf = [sub.wf];
    a = []; for b = 1:length(sub); [m,ii] = max(abs(sub(b).wf));
        if sub(b).wf(ii) < 0; a = [a,b]; end; end
    wf = zscore(wf(:,a));
    shadederrbar(t, nanmean(wf,2), SEM(wf,2), clr{x});
end
xlabel('Time (ms)'); xlim([-2 2]); xticks([-2:2]);
ylabel('Voltage (z-score)'); yticks([-8:4:4]);
legend(lbl,'Location','southeast');
title('Average waveforms'); axis square; set(gca,'TickDir','out');
movegui(gcf,'center');

%% Auto-correlogram
fig = figure; fig.Position([3 4]) = [400 800]; hold on
lbl = {'pSPN','pCIN','other'}; % unit groups
clr = {'k','b','g'}; 
xVal = [-499.85:0.5:-0.35]; xVal = [xVal, -fliplr(xVal)];

for x = unique([units.label])
    idx = find([units.label] == x); % find index for all units of one group
    acgOut = [units(idx).acg]; % extract auto-correlogram
    sm = 2; acgOut = movmean(acgOut,sm,1); % smooth output
    subplot(3,1,x);
    yVal = nanmean(acgOut,2)'; sem = SEM(acgOut,2)';
    color = char2rgb(clr{x}); patchcolor = color+(1-color)*.8;
    yerru = yVal+sem; yerrl = yVal-sem;
    xpatch=[xVal,fliplr(xVal)]; ypatch=[yerru,fliplr(yVal)]; ypatch2=[yVal,zeros(1,length(yVal))];
    hold on
    fill(xpatch,ypatch,patchcolor,'FaceAlpha',0.5,'EdgeAlpha',0);
    fill(xpatch,ypatch2,patchcolor,'FaceAlpha',0.85,'EdgeAlpha',0);
    main = plot(xVal,yVal,'-','Color',color); hold off
    %shadederrbar(xVal, nanmean(acgOut,2), nanstd(acgOut,[],2), 'k');
    xlabel('Lag (ms)'); xlim([-300 300]); xticks([-500:200:500]);
    ylabel('Spikes/s');
    title(lbl{x});
end
movegui(gcf,'center');

%% Comparing properties between units
fig = figure; fig.Position(3) = 700;
lbl = {'pSPN','pCIN','other'}; % unit groups
clr = {'k','b','g'};
lblPlot = {'Firing rate (sp/s)','Coefficient of variation','Phasic activity index','Spike duration (ms)'};
vals = {[units.fr],[units.CV],[units.phasicIndex],[units.spkDur]};
win = [0 50; 0 10; 0 1; 0 4];
bin = [0.5, 0.2, 0.02, 0.1];
for y = 1:4
    plotme = vals{y};
    edges = [win(y,1):bin(y):win(y,2)];
    subplot(2,2,y); hold on
    for x = unique([units.label])
        idx = find([units.label] == x); % find index for all units of one group
        histogram(plotme(idx), edges, 'Normalization', 'probability', 'FaceAlpha', 0.2, 'FaceColor', clr{x}, 'EdgeAlpha', 0.2);
    end
    xlabel(lblPlot{y}); ylabel('% of units');
    title(lblPlot{y}); set(gca,'TickDir','out');
end
movegui(gcf,'center');