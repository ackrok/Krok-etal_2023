load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig1_beh.mat')

%% Extract standard deviation
idxStates = extractBehavioralStates(beh); % identify indices pertaining to different behavioral states
sig = cell(1,2);
for x = 1:length(beh)
    fp = [beh(x).FP{1}, beh(x).FP{2}]; % extract photometry signal
    fp = fp - nanmean(fp);
    for y = 1:2
        sig{y}(x,1) = nanstd(fp(idxStates{x,1},y)); % standard deviation of photometry signal during immobility
        sig{y}(x,2) = nanstd(fp(idxStates{x,2},y)); % standard deviation of photometry signal during locomotion
    end
end
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); nAn = length(uni); % number of unique animals
sig_an = cell(1,2);
for x = 1:nAn
    ii = find(strcmp(rec,uni{x}));
    for y = 1:2
        sig_an{y}(x,:) = nanmean(sig{y}(ii,:),1); % average across multiple recordings from same mouse
    end
end
fprintf('Extracting standard deviation DONE\n');

%%
fig = figure; fig.Position(3) = 1000;
clr = {'g','m'}; lbl = beh(1).FPnames;
for y = 1:2
    subplot(1,2,y); hold on
    plot(sig_an{y}', '.-', 'Color', [0 0 0 0.1], 'MarkerSize', 20);
    errorbar([0.75 2.25], nanmean(sig_an{y}), SEM(sig_an{y},1), '.', 'MarkerSize', 20, 'Color', clr{y});
    xticks([1 2]); xticklabels({'imm','loc'}); 
    ylabel('standard deviation'); ylim([0 6]);
    [~,p] = ttest(sig_an{y}(:,1),sig_an{y}(:,2));
    title(sprintf('%s STD (p = %1.4f) n = %d',lbl{y},p,nAn)); axis square
end
movegui(gcf,'center');