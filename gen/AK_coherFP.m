function [coher_achda, phase_achda, t, f, coher_shuff, phase_shuff] = AK_coherFP(beh, varargin)
%Compute coherogram between ACh and DA photometry signals acquired during
%dual color photometry recordings using Buzsaki function 'MTCoherogram'
%
%   [coher_achda, phase_achda, t, f] = AK_coherFP(beh)
%   [coher_achda, phase_achda, t, f] = AK_coherFP(beh, win)
%   [coher_achda, phase_achda, t, f, coher_shuff, phase_shuff] = AK_coherFP(beh)
%
%   Description: This function is for running cross-correlation analysis on
%   two continuous photometry signals (ACh, DA) after isolating time
%   periods when animal is immobile (no movement, no reward) using the
%   MATLAB function 'xcorr'
%
%   INPUTS
%   'beh' - structure with photometry and behavioral data for multiple
%   recordings, should include beh(x).rec as [an,'-',day]
%   'win'(optional) - window to restrict analysis to, in seconds
%   
%   OUTPUTS
%
%   Author: Anya Krok, March 2022

%% INPUTS
if nargin == 2; win = varargin{1}; end % Window for analysis, in seconds
aa = [2 1]; % Photometry signal to use as reference is the one listed first
% e.g. if y = [1 2] then the signal in beh(x).FP{1} will be used as
% reference signal, while if y = [2 1] then the signal in beh(x).FP{2} will
% be used as reference signal instead
% Default from Krok 2022 is to use y = [2 1] so that rDA1m photometry
% signal is the reference signal
nStates = 3; % Number of behavioral states

%%
mat = struct;
h = waitbar(0,'coherogram');
idxStates = extractBehavioralStates(beh);
for x = 1:length(beh)  % iterate over all recordings
  
    %% extract signals
    mat(x).rec = beh(x).rec; % load recording name
    fp_mat = [];
    fp_mat(:,1) = beh(x).FP{aa(1)}; Fs = beh(x).Fs; % extract photometry signal from structure, which will be used as reference
    fp_mat(:,1) = fp_mat(:,1) - nanmean(fp_mat(:,1)); % subtract baseline (mean of entire photometry signal) from fp
    fp_mat(:,2) = beh(x).FP{aa(2)}; % extract photometry signal from structure
    fp_mat(:,2) = fp_mat(:,2) - nanmean(fp_mat(:,2)); % subtract baseline (mean of entire photometry signal) from fp
    
    %% adjust indices to retain if within specified window
    if nargin == 2
        for z = 1:nStates
            idxTmp = idxStates{x,z};
            idxTmp = idxTmp(idxTmp > win(1)*Fs & idxTmp < win(2)*Fs); % retain only indices that are within specified window
            if ~isempty(idxTmp)
                needL = 200.*Fs;
                if length(idxTmp) < needL % ensure that have at least 200 seconds of signal
                    idxTmp = repmat(idxTmp, [ceil(needL/length(idxTmp)) 1]); % duplicate indices to lengthen signal for processing
                end
            end
            idxStates{x,z} = idxTmp; % re-insert into output structure
        end
    end
    
    %%
    for z = 1:nStates % iterate over behavioral states
        if ~isempty(idxStates{x,z})
            sig = fp_mat(idxStates{x,z},:); % extract indexes samples 
            [coher,ph,t,f] = bz_MTCoherogram(sig(:,1),sig(:,2),'frequency',Fs,'range',[0 10],'window',10,'overlap',5,'step',5,'tapers',[3 5],'pad',0);
            mat(x).coher(:,z) = nanmean(coher,2); % collapse time dimension
            mat(x).phase(:,z) = nanmean(ph,2); % collapse time dimension
        end
    end
    %%
    tmp_coher = []; tmp_phase = []; 
    fp_shuff = fp_mat(:,2);
    for s = 1:50 % repeat shuffle N times
        % fp_shuff = circshift(fp_shuff, Fs);
        fp_shuff = fp_shuff(randperm(length(fp_shuff)));
        [coher,ph] = bz_MTCoherogram(fp_mat(:,1),fp_shuff,'frequency',Fs,'range',[0 10],'window',10,'overlap',5,'step',5,'tapers',[3 5],'pad',0);
        tmp_coher(:,s) = nanmean(coher,2); % collapse time dimension
        tmp_phase(:,s) = nanmean(ph,2); % collapse time dimension
    end
    mat(x).coher_shuff = prctile(tmp_coher, [5 50 95], 2); % 5th, 50th, 95th percentiles
    mat(x).phase_shuff = prctile(tmp_phase, [5 50 95], 2); % 5th, 50th, 95th percentiles
    
    waitbar(x/length(beh),h);
end
close(h); fprintf('Coherogram Analysis Done !\n');

%% EXTRACT COHERENCE and PHASE DURING EACH BEHAVIORAL STATE
coher_state = cell(nStates,1); 
phase_state = cell(nStates,1);
for x = 1:length(mat) % iterate over all recordings
    for z = 1:nStates % iterate over behavioral states
        if z == 3 && size(mat(x).coher,2) < 3 % if no coherence or phase for reward
            coher_state{z}(:,x) = nan(103,1); 
            phase_state{z}(:,x) = nan(103,1);
        else
        coher_state{z}(:,x) = mat(x).coher(:,z); % concatenate
        phase_state{z}(:,x) = mat(x).phase(:,z);
        end
    end
end

%% AVERAGE ACROSS ALL RECORDINGS FOR ONE ANIMAL SUCH THAT N = X mice
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); 
nAn = length(uni); % number of unique animal IDs

coher_achda = cell(3,1); phase_achda = cell(3,1); % initialize output
coher_shuff = cell(3,1); phase_shuff = cell(3,1); % initialize output

for x = 1:nAn
    ii = find(strcmp(rec,uni{x})); % index of matching recordings for this animal
    coh_tmp = [mat(ii).coher_shuff];
    ph_tmp = [mat(ii).phase_shuff];
    for z = 1:nStates
        coher_achda{z}(:,x) = nanmean(coher_state{z}(:,ii),2); % average recordings for each animal, for each behavioral state
        phase_achda{z}(:,x) = nanmean(phase_state{z}(:,ii),2); % average recordings for each animal, for each behavioral state
        coher_shuff{z}(:,x) = nanmean(coh_tmp(:,[z:3:size(coh_tmp,2)]),2); % average recordings for each animal, for each shuff percentile
        phase_shuff{z}(:,x) = nanmean(ph_tmp(:,[z:3:size(ph_tmp,2)]),2); % average recordings for each animal, for each shuff percentile
    end
end
