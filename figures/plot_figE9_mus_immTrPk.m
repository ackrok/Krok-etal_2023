load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_figSX_musc_AChimmTrPk_data.mat');

%% PLOT
if length(cannula) > 2
    choice = menu('Select infusion to compare to saline',{cannula.inf});
    choice = [1, choice];
else
    choice = [1 2];
end

fig = figure; fig.Position([3 4]) = [1375 800];
lbl1 = {'Frequency (events/s)','Duration (ms)','Amplitude (dF/F)'};
lims1 = [0 1.5; 0 700; 0 8];
lbl2 = {'Trough','Peak'};
for z = 1:2 % plot troughs and peaks
    plotme = {[cannula(choice(1)).freq(:,z), cannula(choice(2)).freq(:,z)],...
        [cannula(choice(1)).dur(:,z), cannula(choice(2)).dur(:,z)],...
        [cannula(choice(1)).amp(:,z), cannula(choice(2)).amp(:,z)]};
    for y = 1:length(plotme)
        mat = plotme{y};
        mat = abs(mat);
        % mat(isnan(mat)) = 0;
        % mat(4,:) = nan; % infusion did not work
        % mat = mat./nanmean(mat(:,1));
        [~,p] = ttest(mat(:,1),mat(:,2));
        
        subplot(2,3,y + (z-1)*length(plotme)); hold on
        plot([1 2],mat','.:b','MarkerSize',20);
        errorbar([0.75 2.25],nanmean(mat),SEM(mat,1),'.k','MarkerSize',20);
        xlim([0.5 2.5]); xticks([1 2]); xticklabels({cannula(choice).inf});
        ylabel(sprintf('%s',lbl1{y})); ylim(lims1(y,:));
        title(sprintf('Imm %s %s(p = %1.3f)',lbl2{z},strtok(lbl1{y},'('),p));
        axis('square');
    end
end
movegui(gcf,'center');
