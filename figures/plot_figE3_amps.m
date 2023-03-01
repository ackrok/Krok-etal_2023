load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_fig1_beh.mat')
beh = beh(1:14);

%% INPUTS
choice = menu('','ACh trough + DA peak','ACh peak + DA peak','ACh peak + DA trough');
nFP = 2;
nStates = 3;

%% ANALYSIS
bin = 0.5;
edges = [-15 : bin : 45]; % edges of full spectrum of photometry values
mat = struct;
for x = 1:length(beh)
    %%
    fp_mat = [beh(x).FP{1}, beh(x).FP{2}]; 
    fp_mat = fp_mat - nanmean(fp_mat);
    Fs = beh(x).Fs;
    %
    if isfield(beh,'reward')
        rewYes = extractRewardedTrials(beh(x).reward(:)/Fs, beh(x).lick(:)/Fs, [0 1]);
        rew = beh(x).reward(rewYes);
        idx_rew = extractEventST([1:length(fp_mat)]', floor(rew), floor(rew)+49, 1); % identify sample during reward
    else; idx_rew = []; end
    idx_mov = extractEventST([1:length(fp_mat)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
    idx_imm = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    idx_c = cell(3,1); idx_c{1} = idx_imm_nonRew; idx_c{2} = idx_mov_nonRew; idx_c{3} = idx_rew; % index into cell array for ease of iteration
    %
    binMat = cell(3,1); % chop data into 1s bins
    for a = 1:nStates
        binNum = floor(length(idx_c{a})/Fs);
        sig = fp_mat(idx_c{a},1); % extract ACh signal
        sigBin = reshape(sig(1:Fs*binNum), [Fs, binNum]); % reshape data into 1s bins
        switch choice
            case 1
                sigBin_1 = min(sigBin); % minimum ACh in each 1s bin
            case {2,3}
                sigBin_1 = max(sigBin); % maximum ACh
        end
        sig = fp_mat(idx_c{a},2); % extract DA signal
        sigBin = reshape(sig(1:Fs*binNum), [Fs, binNum]); % reshape data into 1s bins
        switch choice
            case {1,2}
                sigBin_2 = max(sigBin); % maximum DA in each 1s bin
            case 3
                sigBin_2 = min(sigBin); % minimum DA
        end
        binMat{a} = [sigBin_1(:), sigBin_2(:)]; % save ACh minima, DA maxima
    end
    mat(x).rec = beh(x).rec; % add to output structure
    mat(x).binMat = binMat; % add to output structure
end
%
% BY ANIMAL
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); 
nAn = length(uni); % number of unique mice
bin_an = cell(nAn,3); % initialize output
for x = 1:nAn
    ii = find(strcmp(rec,uni{x})); % index of matching recordings for this animal
    for a = 1:nStates
        bin_tmp = [];
        for b = 1:length(ii)
            bin_tmp = [bin_tmp; mat(ii(b)).binMat{a,1}]; % concatenate output for multiple recordings for the same animal
        end
        bin_an{x,a} = bin_tmp;
    end
end
%
% NORMALIZE
bin_an_norm = cell(nAn,nStates);
bin_dist_norm = cell(nAn,nFP);
mode_rew = [];
mid = edges(2:end) - bin/2;
for x = 1:nAn
    for y = 1:nFP
        vals_dist = [];
        for a = 1:nStates
            vals_raw = bin_an{x,a}(:,y); % extract vector
            n = histcounts(vals_raw, edges, 'Normalization', 'probability'); % distribution within specified edges
            vals_dist = [vals_dist, n(:)]; % concatenate output
        end
        [~,ii] = max(vals_dist(:,3)); mode_rew(x,y) = mid(ii); % Mode of distribution for reward
        % Normalize s.t. mode of distribution for reward is 1 and minimum
        % is 0. For ACh, mode is 0 and maximum is 1.
        bin = 0.025;
        edges_norm = [-3 : bin : 3]; % edges of normalized photometry values
        mid_norm = edges_norm(2:end)-bin/2; 
        vals_dist_norm = [];
        for a = 1:nStates
            vals_raw = bin_an{x,a}(:,y);
            vals_norm = vals_raw./abs(mode_rew(x,y)); % Normalize by mode of reward distribution
            bin_an_norm{x,a}(:,y) = vals_norm;
            n = histcounts(vals_norm, edges_norm, 'Normalization', 'probability');
            n = normalize(n, 'range');
            vals_dist_norm = [vals_dist_norm, n(:)];
        end
        bin_dist_norm{x,y} = vals_dist_norm;
    end
end
%
% RE-ORGANIZE OUTPUT
bin_dist_norm_2 = cell(nFP,nStates); 
for y = 1:nFP
    for a = 1:nStates
        for x = 1:nAn
            bin_dist_norm_2{y,a}(:,x) = bin_dist_norm{x,y}(:,a);
        end
    end
end
%
fprintf('Figure S1: amplitudes DONE \n');

%% PLOTTING NORMALIZED
fig = figure; fig.Position(3) = 1000;
clr = {'r','g','b'}; sm = 10;
for y = 1:nFP
    subplot(1,2,y); hold on
    for a = 1:nStates
        shadederrbar(mid_norm, movmean(nanmean(bin_dist_norm_2{y,a},2),sm), movmean(SEM(bin_dist_norm_2{y,a},2),sm), clr{a});
    end
    ylabel('probability'); ylim([-0.05 0.8]); axis('square');
    xlabel('photometry (norm.)');
end
switch choice
    case 1
        subplot(1,2,1); title('ACh minima'); xlim([-3 1]);
        subplot(1,2,2); title('DA maxima'); xlim([-1 3])
    case 2
        subplot(1,2,1); title('ACh maxima'); xlim([-1 3]);
        subplot(1,2,2); title('DA maxima'); xlim([-1 3])
    case 3
        subplot(1,2,2); title('DA minima'); xlim([-3 1]);
        subplot(1,2,1); title('ACh maxima'); xlim([-1 3])
end
movegui(gcf,'center');
% normalized to mode of reward distribution maximuma/minima
% lots of overlap still!

%% PLOT ANIMAL EXAMPLE
x = 11; % example animal

fig = figure; fig.Position([3 4]) = [1375 800]; clearvars sp
colors = [1 0 0 0.1; 0 1 0 0.1; 0 0 1 0.1];
switch choice
    case 1; lbl_ex = {'ACh minima','DA maxima'};
    case 2; lbl_ex = {'ACh maxima','DA maxima'};
    case 3; lbl_ex = {'ACh maxima','DA minima'};
end
bin = 0.5;
edges = [-15 : bin : 45]; % edges of full spectrum of photometry values
mid = edges(2:end) - bin/2;
sm = 5;

for y = 1:2
    subplot(2,3,y); hold on
    for a = 1:3
        n = histcounts(bin_an{x,a}(:,y), edges, 'Normalization', 'probability');
        plot(mid, movmean(n,sm), 'Color', colors(a,[1:3]));
    end
    xlabel(lbl_ex{y}); ylabel('probability'); 
    title(sprintf('%s',lbl_ex{y})); axis square;
end
sp(4) = subplot(2,3,3); hold on
for a = 1:3
    plot(bin_an{x,a}(:,1), bin_an{x,a}(:,2), '.', 'color', colors(a,:), 'MarkerSize', 10);
    title(sprintf('%s',uni{x})); axis equal; axis square;
    xlabel(lbl_ex{1}); ylabel(lbl_ex{2});
end
for a = 1:3
    sp(a) = subplot(2,3,a+3); hold on 
    plot(bin_an{x,a}(:,1), bin_an{x,a}(:,2), '.', 'color', colors(a,:), 'MarkerSize', 10);
    xlabel(lbl_ex{1}); ylabel(lbl_ex{2}); axis square;
end
linkaxes(sp,'x'); linkaxes(sp,'y'); 
switch choice
    case 1
        xlim([-8 2]); ylim([-5 17]);
        subplot(2,3,1); xlim([-8 2]); subplot(2,3,2); xlim([-5 17]);
    case 3
        xlim([-5 30]); ylim([-8 4]);
        subplot(2,3,1); xlim([-5 30]); subplot(2,3,2); xlim([-8 4]); 
end
movegui(gcf,'center');