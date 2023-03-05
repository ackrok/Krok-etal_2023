load('C:\Users\Anya\Desktop\FP_LOCAL\krok_fig1_beh.mat')

%% Figure 1B
flbl = 'fig 1B'; state = {'reward','locomotion','immobility'};
a = [20 20 17]; % for plotting
c = [751.5 755; 249.5 253; 1559.5 1563]; % x-lim for plotting
fig = figure; fig.Position(3) = 1375;
for b = 1:3
    sp(b) = subplot(1,3,b); hold on; 
    x = a(b);
    fp_mat = [beh(x).FP{1}, beh(x).FP{2}]; % extract photometry signal from structure
    fp_mat = fp_mat - nanmean(fp_mat); % subtract mean, now centered on zero
    Fs = beh(x).Fs; % sampling frequency
    plot(beh(x).time, fp_mat(:,1), 'g'); % plot ACh photometry signal
    plot(beh(x).time, fp_mat(:,2), 'm'); % plot DA photometry signal
    if ~isempty(beh(x).reward)
        stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1),'b'); % plot reward delivery times
        [~,~,lick] = extractRewardedTrials(beh(x).reward./Fs, beh(x).lick./Fs);
        plot([lick,lick]',[10;11].*ones(2,length(lick)),'-k'); % plot lick times
    end
    plot(beh(x).time, getAcc(beh(x).vel) - 5, 'k'); % plot acceleration
    xlabel('time (s)'); xlim([c(b,:)]);
    ylabel('fluorescence (%dF/F)');
    title(sprintf('%s - %s',flbl,state{b}));
    axis square
end
linkaxes(sp,'y');

%% Figure S1C
flbl = 'fig S1C'; state = {'reward','locomotion','immobility'};
a = [23]; % for plotting
c = [1573 1596]; % x-lim for plotting
fig = figure; fig.Position(3) = 1375;
for b = 1
    hold on; 
    x = a(b);
    fp_mat = [beh(x).FP{1}, beh(x).FP{2}]; % extract photometry signal from structure
    fp_mat = fp_mat - nanmean(fp_mat); % subtract mean, now centered on zero
    Fs = beh(x).Fs; % sampling frequency
    plot(beh(x).time, fp_mat(:,1), 'g'); % plot ACh photometry signal
    plot(beh(x).time, fp_mat(:,2), 'm'); % plot DA photometry signal
    if ~isempty(beh(x).reward)
        stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1),'b'); % plot reward delivery times
        [~,~,lick] = extractRewardedTrials(beh(x).reward./Fs, beh(x).lick./Fs);
        plot([lick,lick]',[10;11].*ones(2,length(lick)),'-k'); % plot lick times
    end
    plot(beh(x).time, getAcc(beh(x).vel) - 5, 'k'); % plot acceleration
    xlabel('time (s)'); xlim([c(b,:)]);
    ylabel('fluorescence (%dF/F)');
    title(sprintf('%s - %s',flbl,state{b}));
end

%% Figure S1D
flbl = 'fig S1D'; state = {'reward','locomotion','immobility'};
a = [21 18 18]; % for plotting
c = [939 944; 22 27; 1091.5 1096.5]; % x-lim for plotting
fig = figure; fig.Position(3) = 1375;
for b = 1:3
    sp(b) = subplot(1,3,b); hold on; 
    x = a(b);
    fp_mat = [beh(x).FP{1}, beh(x).FP{2}]; % extract photometry signal from structure
    fp_mat = fp_mat - nanmean(fp_mat); % subtract mean, now centered on zero
    Fs = beh(x).Fs; % sampling frequency
    plot(beh(x).time, fp_mat(:,1), 'g'); % plot ACh photometry signal
    plot(beh(x).time, fp_mat(:,2), 'm'); % plot DA photometry signal
    if ~isempty(beh(x).reward)
        stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1),'b'); % plot reward delivery times
        [~,~,lick] = extractRewardedTrials(beh(x).reward./Fs, beh(x).lick./Fs);
        plot([lick,lick]',[10;11].*ones(2,length(lick)),'-k'); % plot lick times
    end
    plot(beh(x).time, getAcc(beh(x).vel) - 5, 'k'); % plot acceleration
    xlabel('time (s)'); xlim([c(b,:)]);
    ylabel('fluorescence (%dF/F)');
    title(sprintf('%s - %s',flbl,state{b}));
    axis square
end
linkaxes(sp,'y');

%% Figure S1E
flbl = 'fig S1E'; state = {'reward','locomotion','immobility'};
a = [15 15 15]; % for plotting
c = [280 285; 1594 1599; 1671 1676]; % x-lim for plotting
fig = figure; fig.Position(3) = 1375;
for b = 1:3
    sp(b) = subplot(1,3,b); hold on; 
    x = a(b);
    fp_mat = [beh(x).FP{1}, beh(x).FP{2}]; % extract photometry signal from structure
    fp_mat = fp_mat - nanmean(fp_mat); % subtract mean, now centered on zero
    Fs = beh(x).Fs; % sampling frequency
    plot(beh(x).time, fp_mat(:,1), 'g'); % plot ACh photometry signal
    plot(beh(x).time, fp_mat(:,2), 'm'); % plot DA photometry signal
    if ~isempty(beh(x).reward)
        stem(beh(x).reward./Fs, 10.*ones(length(beh(x).reward),1),'b'); % plot reward delivery times
        [~,~,lick] = extractRewardedTrials(beh(x).reward./Fs, beh(x).lick./Fs);
        plot([lick,lick]',[10;11].*ones(2,length(lick)),'-k'); % plot lick times
    end
    plot(beh(x).time, getAcc(beh(x).vel) - 5, 'k'); % plot acceleration
    xlabel('time (s)'); xlim([c(b,:)]);
    ylabel('fluorescence (%dF/F)');
    title(sprintf('%s - %s',flbl,state{b}));
    axis square
end
linkaxes(sp,'y');

