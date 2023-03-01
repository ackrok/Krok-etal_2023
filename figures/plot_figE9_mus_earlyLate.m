load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_snc_salVmus_t0v700s_data');

%%
r = [0.02; 700]; r = [r, r+10];
val = cell(2,2);
for z = 1:2
    beh = cannula(z).s; beh(1) = []; nAn = length(beh);
    for x = 1:length(beh)
        Fs = beh(x).Fs;
        for y = 1:2
            for ii = 1:2
                val{z,y}(x,ii) = nanmean(beh(x).nbFP{y}(r(ii,1)*Fs : r(ii,2)*Fs));
            end
        end
    end
end

fig = figure; fig.Position(3) = 1000;
clr = {'g','m'};
for y = 1:2
    subplot(1,2,y); hold on
    plot([1;2].*ones(2,nAn), val{1,y}','.-k','MarkerSize',20);
    plot([3;4].*ones(2,nAn), val{2,y}','.-','MarkerSize',20,'Color',clr{y});
    xlim([0.5 4.5]); xticks([1:4]); xticklabels({'SAL 0','SAL X','MUS 0','MUS X'});
    ylabel('a.u.');
    p = []; for z = 1:2; [~,p(z)] = ttest(val{z,y}(:,1),val{z,y}(:,2)); end
    title(sprintf('t0 vs %d min | SAL (p = %1.4f), MUS (p = %1.3f)',round(r(2,1)/60),p(1),p(2)));
end

%%
r = [0.02; 700]; r = [r, r+10];
val = cell(2,2);
for z = 1:2
    beh = cannula(z).s; beh(1) = []; nAn = length(beh);
    for x = 1:length(beh)
        Fs = beh(x).Fs;
        for y = 1:2
            for ii = 1:2
                val{z,y}(x,ii) = nanmean(beh(x).nbFP{y}(r(ii,1)*Fs : r(ii,2)*Fs));
            end
        end
    end
end


fig = figure; fig.Position(3) = 1000;
clr = {'g','m'};
z = 0; j1 = z-0.3; j2 = z+0.3; jit = j1 + (j2-j1).*rand(nAn,1);
for y = 1:2
    subplot(1,2,y); hold on
    z = 1; plot(jit+z, val{z,y}(:,2)./val{z,y}(:,1), '.','MarkerSize',20,'Color','k');
    z = 2; plot(jit+z, val{z,y}(:,2)./val{z,y}(:,1), '.','MarkerSize',20,'Color',clr{y});
    a = [val{1,y}(:,2)./val{1,y}(:,1), val{2,y}(:,2)./val{2,y}(:,1)];
    errorbar([1 2], nanmean(a), SEM(a,1), '.-r', 'MarkerSize',20);
    xlim([0.5 2.5]); xticks([1:2]); xticklabels({'SAL','MUS'});
    ylabel('a.u. (prop. relative to t0)'); ylim([0.5 1.1])
    [~,p2] = ttest(a(:,1),a(:,2),'tail','right');
    title(sprintf('%d min | SAL vs MUS (p = %1.4f)',round(r(2,1)/60),p2));
end

%%
bb = [100:50:1500];
p2 = [];
for aa = 1:length(bb)
r = [0.02; bb(aa)]; r = [r, r+10];
val = cell(2,2);
for z = 1:2
    beh = cannula(z).s; beh(1) = []; nAn = length(beh);
    for x = 1:length(beh)
        Fs = beh(x).Fs;
        for y = 1:2
            for ii = 1:2
                val{z,y}(x,ii) = nanmean(beh(x).nbFP{y}(r(ii,1)*Fs : r(ii,2)*Fs));
            end
        end
    end
end

y = 1;
a = [val{1,y}(:,2)./val{1,y}(:,1), val{2,y}(:,2)./val{2,y}(:,1)];
[~,p2(aa)] = ttest(a(:,1),a(:,2));
end

figure;
plot(bb, p2, '.k'); 
ylabel('p-value SAL vs MUS');
xlabel('tX to compare to t0')
title('ACh t0 vs tX');

%%
bb = [100:50:1500];
p2 = [];
for aa = 1:length(bb)
r = [0.02; bb(aa)]; r = [r, r+10];
val = cell(2,2);
for z = 1:2
    beh = cannula(z).s; beh(1) = []; nAn = length(beh);
    for x = 1:length(beh)
        Fs = beh(x).Fs;
        for y = 1:2
            for ii = 1:2
                val{z,y}(x,ii) = nanmean(beh(x).nbFP{y}(r(ii,1)*Fs : r(ii,2)*Fs));
            end
        end
    end
end

y = 2;
[~,p2(aa)] = ttest(val{1,y}(:,1),val{1,y}(:,2)); 
end

figure;
plot(bb, p2, '.k'); 
ylabel('p-value SAL t0 vs SAL tX');
xlabel('tX to compare to t0 (s)')
title('DA - SAL t0 vs tX');