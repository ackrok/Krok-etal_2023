%%
load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_SX_globalSignal_data');

%%
r = [0.02; 1200]; r = [r, r+100];
val = cell(1,2); for y = 1:2; val{y} = nan(8,length(cannula)); end
mouse = {'AK189','AK190','AK193','AK197','AK239','AK240','AK241','AK243'};

for z = 1:length(cannula)
    beh = cannula(z).s; nAn = length(beh);
    for x = 1:length(beh)
        Fs = beh(x).Fs;
        idx = strcmp(mouse,strtok(beh(x).rec,'-'));
        for y = 1:2
            tmp = [];
            for ii = 1:2
                tmp(ii) = nanmean(beh(x).nbFP{y}(r(ii,1)*Fs : r(ii,2)*Fs));
            end
            val{y}(idx,z) = tmp(2)/tmp(1);
        end
    end
end

%
fig = figure; fig.Position(3) = 1000;
clr = {'k','r','g','m','b'}; lbl = {cannula.inf};
z = 0; j1 = z-0.3; j2 = z+0.3; jit = j1 + (j2-j1).*rand(nAn,1);
for y = 1:2
    subplot(1,2,y); hold on
    for z = 1:length(cannula)
        plot(jit+z, val{y}(:,z), '.','MarkerSize',20,'Color',clr{z}); end
    % a = [val{1,y}(:,2)./val{1,y}(:,1), val{2,y}(:,2)./val{2,y}(:,1)];
    % errorbar([1 2], nanmean(a), SEM(a,1), '.-r', 'MarkerSize',20);
    xlim([0.5 5.5]); xticks([1:5]); xticklabels(lbl);
    ylabel('a.u. (prop. relative to t0)'); ylim([0 1.1])
    % [~,p2] = ttest(a(:,1),a(:,2),'tail','right');
    % title(sprintf('%d min | SAL vs MUS (p = %1.4f)',round(r(2,1)/60),p2));
end

%%
figure; 
for x = 1:4
    subplot(2,2,x); hold on; 
    clr = {'k','r','g','m'};
    for z = [1 2 4]
        plot(cannula(z).s(1).time, cannula(z).s(1).nbFP{2},clr{z}); 
    end
    legend({'saline DA','iGluR DA','D1/2R DA'}); 
end
%%
figure; 
for x = 5:8
    subplot(2,2,x-4); hold on; 
    clr = {'k','r','g','m'};
    for z = [1 2 3]
        plot(cannula(z).s(1).time, cannula(z).s(1).nbFP{1},clr{z}); 
    end
    legend({'saline ACh','iGluR ACh','mAChR ACh'}); 
end
