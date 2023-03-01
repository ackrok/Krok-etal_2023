%% LOAD
load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig4_st2ach');
% IV042 rec05, IV043 rec01, IV043 rec02, IV046 rec02, IV046 rec04, IV066
% rec01, IV066 rec03

%% (2) PLOT
for x = 1:2 % firing rate to ACh peaks and troughs
    fig = figure; fig.Position(3) = 1375;
    time = out(x).time;
    nUnits = size(out(x).delta,2);
    lims_1 = [3 10; 2 8]; lims_2 = [0.6 1.8; 0.6 1.3]; lims_3 = [0.8 1.4; 0.7 1.1];

    subplot(1,3,1); hold on
    sm = 5;
    plot(out(x).time, movmean(out(x).example,5,1));
    plot([0 0],lims_1(x,:),'--','Color',[0.5 0.5 0.5]);
    ylabel('pCIN firing rate (sp/s)'); ylim(lims_1(x,:));
    xlabel(sprintf('Time to %s (s)',out(x).lbl)); xlim([-1 1]); xticks([-1:0.5:1]);
    title('pCIN example (14, 16) immobility');
    axis('square'); set(gca,'TickDir','out');

    subplot(1,3,2); hold on
    switch x
        case 1; [~, ii] = sort(max(out(x).delta)); % sort in ascending order
        case 2; [~, ii] = sort(min(out(x).delta)); % sort in ascending order
    end
    h = imagesc(time, [1:nUnits], 1+out(x).delta(:,ii)', lims_2(x,:));
    colorbar; colormap(jet(256));
    xlabel(sprintf('Time to %s (s)',out(x).lbl)); xlim([-1 1]); xticks([-1:0.5:1]);
    ylabel('Unit');
    h(x) = h.Parent; h(x).CLim = lims_2(x,:);
    title(sprintf('immobility (n = %d)',nUnits));
    axis square; set(gca,'TickDir','out');

    subplot(1,3,3); hold on
    delta = out(x).delta - nanmean(out(x).delta([1:find(time == -0.51)],:)); % subtract baseline
    shuff = out(x).delta50 - nanmean(out(x).delta50([1:find(time == -0.51)],:)); % subtract baseline
    plot([0 0],lims_3(x,:),'--','Color',[0.5 0.5 0.5]);
    shadederrbar(time, 1+movmean(nanmean(shuff,2),sm), movmean(nanmean(out(x).delta95,2),sm), 'k');
    shadederrbar(time, 1+movmean(nanmean(delta,2),1), movmean(SEM(delta,2),1), 'b'); 
    ylabel('Firing rate (norm.)'); ylim(lims_3(x,:)); yticks([0:0.1:2]);
    xlabel(sprintf('Time to %s (s)',out(x).lbl)); xlim([-1 1]); xticks([-1:0.5:1]);
    switch x
        case 1; title(sprintf('max = %1.3f',max(1+nanmean(out(x).delta,2))));
        case 2; title(sprintf('min = %1.3f',min(1+nanmean(out(x).delta,2))));
    end
    axis('square'); set(gca,'TickDir','out');

    movegui(gcf,'center');

    % subplot(1,3,3); hold on
    % ds = 2;
    % a = 100*sum(out(x).above95,2)/size(out(x).above95,2); 
    % [prop_m,prop_t] = max(a); prop_t = time(prop_t); 
    % a = a(2:ds:end);
    % b = -100*sum(out(x).below5,2)/size(out(x).below5,2); b = b(1:ds:end);
    % bar(time(2:ds:end), a,'FaceColor','b','FaceAlpha',0.5,'EdgeAlpha',0.1);
    % bar(time(2:ds:end), b,'FaceColor','b','FaceAlpha',0.5,'EdgeAlpha',0.1);
    % xlabel(sprintf('Time to %s (s)',out(x).lbl)); xlim([-1 1]); xticks([-1:0.5:1]);
    % ylabel('% of units'); ylim([-100 100]);
    % title(sprintf('max %1.1f @ %d ms',prop_m, round(1000*prop_t)))
    % axis('square'); set(gca,'TickDir','out');
    % movegui(gcf,'center');
end