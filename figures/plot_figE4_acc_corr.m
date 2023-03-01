load('krok_fig1_reward_beh.mat')
beh1 = beh;

%%
beh = beh1;
for x = 1:length(beh)
    acc = getAcc(beh1(x).vel);
    % beh(x).FP{2} = acc; % ACh vs Acceleration, with acceleration as reference
    beh(x).FP{1} = beh(x).FP{2}; beh(x).FP{1} = acc; % DA vs Acc
end
[corr_da, lags, shuff_da] = AK_corrFP(beh);

%%
figure; hold on
Fs = beh(1).Fs;
z = 2; % Locomotion only
shadederrbar(lags/Fs, nanmean(corr_ach{z},2), SEM(corr_ach{z},2), 'g');
shadederrbar(lags/Fs, nanmean(corr_da{z},2), SEM(corr_da{z},2), 'm');
plot([0 0],[-0.5 0.5],'--k');
ylim([-0.2 0.4]); ylabel('coefficient');
xlim([-2 2]); xlabel('Lag to acceleration (s)');
legend({'ACh vs acc','DA vs acc'});
title('avg corr, loc only, acc reference'); axis square
