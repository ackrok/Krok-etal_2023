function [fp2ph, fp2ph_norm, da2achph_norm, mid, da2achph] = AK_fp2phase(beh)
%Distribution of signal fluorescence to the instantaneous phase of a
%second photometry signal, acquired using dual color fiber photometry
%
%   [fp2ph, fp2ph_norm, da2achph_norm, mid, da2achph] = AK_fp2phase(beh)
%
%   Description: This function is for computing the distribution of
%   photometry signal to instantaneous phase of the same or another
%   simultaneously recording photometry signal.
%   First, will compute the instantaneous phase of each signal
%   Then for each signal, will determine the fluorescence (%dF/F, normalized)
%   at each instantaneous phase, across behavioral states
%   Lastly, will compute the normalized fluorescence for rDA1m signal at
%   different instantaneous phase of ACh3.0 signal, across behavioral states
%
%   INPUTS
%   'beh' - structure with photometry and behavioral data for multiple
%   recordings, should include beh(x).rec as [an,'-',day]
%   
%   OUTPUTS
%   'fp2ph' - Fluorescence (%dF/F) aligned to instantaneous phase
%       fp2ph{a,b}, where a designates the fluorescence signal (e.g. rDA1m)
%       and b designates the behavioral state
%           a = 1 is for ACh3.0, a = 2 is for rDA1m
%           b = 1 is immobility, b = 2 is locomotion, b = 3 is reward
%   'fp2ph_norm' - Fluorescence (normalized) aligned to instantaneous phase
%       fp2ph{a,b}, same as above for 'fp2ph'
%   'da2achph_norm' - rDA1m normalized fluorescence aligned to
%   instantaneous phase of ACh3.0 fluorescence signal
%       da2achph_norm{a,b}, where a designates which fluorescence signal's
%       instantaneous phase is being aligned to
%       b designates the behavioral state
%           a = 1 is for rDA1m normalized fluorescence aligned to
%           instantaneous phase of ACh3.0 fluorescence
%           a = 2 is for ACh3.0 normalized fluorescence aligned to
%           instantaneous phase of rDA1m fluorescence
%           b = 1 is immobility, b = 2 is locomotion, b = 3 is reward
%   'mid' - Vector of midpoints of bins (-180 to +180 degrees), for plotting
%
%   Author: Anya Krok, March 2022

%% INPUTS
NumStd = 1.5; % Phase = 0 ("peaks") are identified as points that cross NumStd standard deviations (default 1.5)
aa = [1 2]; % Photometry signal to use as reference is the one listed first
% e.g. if y = [1 2] then the signal in beh(x).FP{1} will be used as
% reference signal, while if y = [2 1] then the signal in beh(x).FP{2} will
% be used as reference signal instead
% Default from Krok 2022 is to use y = [1 2] so that ACh3.0 photometry
% signal is the reference signal and rDA1m photometry is thus aligned to
% the instantaneous phase of the ACh3.0 signal 
nStates = 3; % Number of behavioral states

%%
fp2ph_beh = cell(2,3); % initialize termporary output, aligning photometry signal to instananeous phase for each recording
fp2ph_norm_beh = fp2ph_beh;
da2achph_beh = fp2ph_beh;
da2achph_norm_beh = fp2ph_beh;

h = waitbar(0, 'photometry to instantaneous phase');
for x = 1:length(beh)
    Fs = beh(x).Fs;
    fp_mat = []; fp_phase = []; fp_deg = [];
    fp2ph_tmp = cell(2,2); for a = 1:2; for b = 1:2; fp2ph_tmp{a,b} = nan(36, 10000); end; end % initialize matrix
    fp2ph_state_tmp = cell(2,3); for a = 1:2; for b = 1:3; fp2ph_state_tmp{a,b} = nan(36, 10000); end; end % initialize matrix
    fp2ph_state_norm_tmp = fp2ph_state_tmp; % initialize matrix

    for y = 1:2
        fp_mat(:,y) = beh(x).FP{aa(y)}; % extract photometry signal from structure, which will be used as reference
        fp_mat(:,y) = fp_mat(:,y) - nanmean(fp_mat(:,y)); % subtract baseline (mean of entire photometry signal) from fp
        
        %% Index of behavioral states
        if isfield(beh,'reward')
            idx_rew = extractEventST([1:length(fp_mat)]', floor(beh(x).reward), floor(beh(x).reward)+(Fs*2), 1); % identify sample during reward
        else; idx_rew = []; end
        idx_mov = extractEventST([1:length(fp_mat)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
        idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
        idx_imm = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
        idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
        idx_c = cell(nStates,1); idx_c{1} = idx_imm_nonRew; idx_c{2} = idx_mov_nonRew; idx_c{3} = idx_rew; % index into cell array for ease of iteration
        
        %% Bandpass filter
        signal = fp_mat(:,y);
        Fpass = [0.5 4];
        Fs = beh(x).Fs; %sampling rate, has to be at least double of your high pass frequency
        Wn = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; 
        [b,a] = butter(3,Wn);
        data_filt= filtfilt(b,a,signal); % signal is your lfp data, output is the filtered data

        %% Instantaneous phase
        H = hilbert(double(data_filt));
        data_phase = angle(H); % output is the instantaneous phase
        fp_phase(:,y) = data_phase;
        fp_deg(:,y) = rad2deg(data_phase);

        %Take the absolute of the filtered signal and calculate the standard deviation
        rmssig  = abs(data_filt);
        stdsig = std(rmssig);

        %Find the index of peaks that are more than NumStd standard deviations
        % finds the 0 degree phase indices that cross 1.5 standard deviations.
        peakIdx = find(data_filt>(NumStd*stdsig) & [0; diff(data_phase>0)]);         

        %% Align photometry to phase
        bin = 10;
        [~, edges] = histcounts(fp_deg(:,y), 'BinWidth', bin, 'BinLimits', [-180 180]); % genereate edges vector for phase, in degrees
        mid = edges(2:end) - bin/2;
        [phIdx] = discretize(fp_deg(:,y), edges); % returns the indices of the bins that the elements of X fall into

        for z = 1:length(unique(phIdx))
            ii = find(phIdx == z); % find index in signal that are this phase
            fp2ph_tmp{y,1}(z,[1:length(ii)]) = fp_mat(ii,y); % extract photometry %dF/F values matching this phase
            for c = 1:length(idx_c)
                ii_c = ii(ismember(ii, idx_c{c})); % include index for this behavioral state
                fp2ph_state_tmp{y,c}(z,[1:length(ii_c)]) = fp_mat(ii_c,y); % extract photometry values for this behavioral state
            end
        end
        [~,b] = find(isnan(fp2ph_tmp{y,1})); % find first column with NaN's
        fp2ph_tmp{y,1}(:,[b(1):size(fp2ph_tmp{y,1},2)]) = []; % remove end of matrix with incomplete columns
        fp2ph_tmp{y,2} = normalize(fp2ph_tmp{y,1},1,'range'); % normalize over each individual [0 180] oscillation
        % [~,b] = find(isnan(fp2ph{y,2})); % find first column with NaN's
        % fp2ph{y,2}(:,[b(1):size(fp2ph{y,2},2)]) = []; % remove end of matrix with incomplete columns

        for z = 1:length(idx_c)
            fp2ph_state_norm_tmp{y,z} = normalize(fp2ph_state_tmp{y,z},1,'range');
            [~,b] = find(isnan(fp2ph_state_tmp{y,z})); % find first column with NaN's
            fp2ph_state_norm_tmp{y,z}(:,[b(1):size(fp2ph_state_tmp{y,z},2)]) = []; % remove end of matrix with incomplete columns
        end

        %% save into output cell array
        for z = 1:length(idx_c)
            fp2ph_beh{y,z}(:,x) = nanmean(fp2ph_state_tmp{y,z},2);
            fp2ph_norm_beh{y,z}(:,x) = nanmean(fp2ph_state_norm_tmp{y,z},2);
        end
    end

    %% DA fp norm to ACh phase
    for xx = 1:2
        switch xx 
            case 1
                y = 1;
                [~, edges] = histcounts(fp_deg(:,y), 'BinWidth', 10, 'BinLimits', [-180 180]); % genereate edges vector for phase, in degrees
                [phIdx] = discretize(fp_deg(:,y), edges); % returns the indices of the bins that the elements of X fall into
                y = 2;  
            case 2
                y = 2;
                [~, edges] = histcounts(fp_deg(:,y), 'BinWidth', 10, 'BinLimits', [-180 180]); % genereate edges vector for phase, in degrees
                [phIdx] = discretize(fp_deg(:,y), edges); % returns the indices of the bins that the elements of X fall into
                y = 1;      
        end
        da2achph_state = cell(3,1); for a = 1:3; da2achph_state{a} = nan(36, 10000); end % initialize matrix
        da2achph_state_norm = cell(3,1); for a = 1:3; da2achph_state_norm{a} = nan(36, 10000); end % initialize matrix
        for z = 1:length(unique(phIdx))
            ii = find(phIdx == z); % find index in signal that are this phase
            for c = 1:length(idx_c)
                ii_c = ii(ismember(ii, idx_c{c})); % include index for this behavioral state
                da2achph_state{c}(z,[1:length(ii_c)]) = fp_mat(ii_c,y); % extract photometry values for this behavioral state
            end
        end
        for z = 1:length(idx_c)
            da2achph_state_norm{z} = normalize(da2achph_state{z},1,'range');
            [~,b] = find(isnan(da2achph_state_norm{z})); % find first column with NaN's
            da2achph_state_norm{z}(:,[b(1):size(da2achph_state{z},2)]) = []; % remove end of matrix with incomplete columns
            da2achph_norm_beh{xx,z}(:,x) = nanmean(da2achph_state_norm{z},2); % save into output cell array
            da2achph_beh{xx,z}(:,x) = nanmean(da2achph_state{z},2);
        end
    end

    %%
    waitbar(x/length(beh),h);
end
fprintf('Done! \n'); close(h);

%% AVERAGE ACROSS ALL RECORDINGS FOR ONE ANIMAL SUCH THAT N = X mice
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); 
nAn = length(uni); % number of unique animal IDs

fp2ph = cell(2,3); % Initialize cell arrays
fp2ph_norm = fp2ph;
da2achph_norm = fp2ph; 
da2achph = fp2ph; 

for x = 1:nAn
    ii = find(strcmp(rec,uni{x}));
    for z = 1:length(idx_c)
        for y = 1:2
            fp2ph{y,z}(:,x) = nanmean([fp2ph_beh{y,z}(:,ii)],2);
            fp2ph_norm{y,z}(:,x) = nanmean([fp2ph_norm_beh{y,z}(:,ii)],2);
            fp2ph_norm{y,z}(:,x) = normalize(fp2ph_norm{y,z}(:,x), 'range');
            da2achph_norm{y,z}(:,x) = nanmean([da2achph_norm_beh{y,z}(:,ii)],2);
            da2achph_norm{y,z}(:,x) = normalize(da2achph_norm{y,z}(:,x), 'range');
            da2achph{y,z}(:,x) = nanmean([da2achph_beh{y,z}(:,ii)],2);
        end
    end
end
