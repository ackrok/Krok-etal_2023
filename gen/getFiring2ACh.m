function out = getFiring2ACh(units, beh)
% Align firing rate of units to peaks and troughs of ACh photometry signal
%
% out = getFiring2ACh(units, beh)
%
% INPUTS
%   'units' - structure containing spike times for units and recording ID, 
%       for matching to behavioral data in 'beh'
%   'beh' - structure containing photometry signal and behavioral data, 
%       including on and off times for periods of immobility and locomotion
%       and must also include recording ID in matching format for matching 
%       behavioral data to spiking data
%
% OUTPUT
%   'out' - structure containing firing rates aligns to peaks and troughs
%   of photometry signal
%
% Anya Krok, January 2022

out = struct;
lbl = {'peak','trough'};
for zz = 1:2
    out(zz).lbl = lbl{zz};

nShuff = 5;
bin = 0.02; window = [-1 1]; %CHANGE: window for PETH
mat = struct; % Initialize temporary output structure
h = waitbar(0, 'pCIN spikes to ACh peaks/troughs');
for x = 1:length(beh)
    %% Extract spike times
    idx = find(strcmp({units.rec},beh(x).rec));
    if isempty(idx); continue; end
    st = {units(idx).st}; % Extract spike times of units from this recording
    fr = []; % Compute unit firing rate during movement and rest
    Fs = beh(x).Fs; 
    if all(logical(~rem(beh(x).onRest,1))); diffFs = Fs; else; diffFs = 1; end % Not divisible if seconds, then do not adjust values
    for y = 1:length(idx)
        % fr(y,1) = 1/mean(diff(extractEventST(st{y},beh(x).on./diffFs,beh(x).off./diffFs,0))); % Event times during locomotion
        fr(y,2) = 1/mean(diff(extractEventST(st{y},beh(x).onRest./diffFs,beh(x).offRest./diffFs,0))); % Event times during immobility
    end

    %% Identify peaks and troughs of ACh photometry signal
    % Bandpass filter
    signal = beh(x).FP{1} - nanmean(beh(x).FP{1}); % PHOTOMETRY
    Fpass = [0.5 4];
    Fs = 50; %sampling rate, has to be at least double of your high pass frequency
    Wn = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; 
    [b,a] = butter(3,Wn);
    dataFilt = filtfilt(b,a,signal); % signal is your lfp data, output is the filtered data
    H = hilbert(double(dataFilt));
    dataPhase = angle(H); % output is the instantaneous phase
    % Take the absolute of the filtered signal and calculate the standard deviation
    rmssig  = abs(dataFilt);
    stdsig = std(rmssig);
    % Find the index of peaks that are more than NumStd standard deviations (NumStd could be 1.5 standard deviations)
    % finds the 0 degree phase indices that cross 1.5 standard deviations. Might have to remove the ; depending on whether your data is a row or column array
    NumStd = 1.5;
    idxPeak = find(dataFilt>(NumStd*stdsig) & [0; diff(dataPhase>0)]);         
    idxTrough = find(-dataFilt>(NumStd*stdsig) & [0; diff(-dataPhase>0)]);
    switch zz
        case 1; locs = idxPeak./Fs; % convert to seconds
        case 2; locs = idxTrough./Fs; 
    end
    
    %% Extract event times during rest & movement
    evSub = cell(1,2); 
    % evSub{1} = extractEventST(locs, beh(x).on./Fs, beh(x).off./Fs, 1); % Event times during locomotion
    evSub{2} = extractEventST(locs, beh(x).onRest./Fs, beh(x).offRest./Fs, 1); % Event times during rest

        %% PETH
    align = cell(1,length(evSub)); 
    alignDelta = cell(1,length(evSub)); 
    shuff95 = align; shuff50 = align; shuff5 = align;
    for z = 1:length(evSub)
        if isempty(evSub{z}); continue; end
        peth = getClusterPETH(st, evSub{z}, bin, window); % PETH: spike times aligned to fp peaks
        align{z} = peth.fr;
        h2 = waitbar(0, 'Aligning spike times to ACh peaks/troughs');
        for y = 1:length(idx)
            alignDelta{z} = [alignDelta{z}, (peth.fr(:,y)-fr(y,z))./fr(y,z)]; % Delta firing rate change
            stShuff = shuffleST(st{y}, nShuff);
            peth_shuff = getClusterPETH(stShuff, evSub{z}, bin, window); %PETH: shuffled spike times aligned to fp peaks
            prc = prctile(peth_shuff.fr,[2.5 50 97.5],2); %5th, 50th, 95th percentile of shuffled PETH
            shuff50{z} = [shuff50{z}, (prc(:,2)-fr(y,z))./fr(y,z)]; 
            shuff95{z} = [shuff95{z}, (prc(:,3)-fr(y,z))./fr(y,z)]; 
            shuff5{z} = [shuff5{z}, (prc(:,1)-fr(y,z))./fr(y,z)]; 
            waitbar(y/length(idx),h2);
        end
        close(h2);
    end
    
    %% STA
%     sta_fp = cell(1,length(evSub));
%     for z = 1:length(evSub)
%         if isempty(evSub{z}); continue; end
%         [sta_fp{z},sta_time] = getSTA(signal, evSub{z}, Fs, [window(1), window(2)]); % Align photometry to FP peaks
%     end
    
    %% Load into output structure
    mat(x).rec = beh(x).rec; 
%     mat(x).FPnames{1} = beh(x).FPnames{1}; 
%     mat(x).n = [units(idx).n]; 
%     mat(x).fr = fr;
    mat(x).align_rest = align{2}; 
    mat(x).alignDelta_rest = alignDelta{2}; 
    mat(x).shuff5_rest = shuff5{2}; 
    mat(x).shuff50_rest = shuff50{2}; 
    mat(x).shuff95_rest = shuff95{2}; 
%     mat(x).sta_rest = sta_fp{2}; 
   %%
    waitbar(x/length(beh),h); 
end
close(h);
time = peth.time; 

%% Load into output structure
align_rest = []; alignDelta_rest = []; 
shuff50_rest = []; shuff95_rest = []; shuff5_rest = [];
for x = 1:length(mat)
    if isempty(mat(x).align_rest); continue; end
    align_rest = [align_rest, mat(x).align_rest];
    alignDelta_rest = [alignDelta_rest, mat(x).alignDelta_rest];
    shuff95_rest = [shuff95_rest, mat(x).shuff95_rest];
    shuff50_rest = [shuff50_rest, mat(x).shuff50_rest];
    shuff5_rest = [shuff5_rest, mat(x).shuff5_rest];
end
above95 = []; below5 = [];
for y = 1:size(alignDelta_rest,2)
    a = alignDelta_rest(:,y);
    a = a - nanmean(a([1:find(time == -0.51)],:));
    b = shuff95_rest(:,y);
    c = shuff50_rest(:,y) - (shuff95_rest(:,y) - shuff50_rest(:,y));  % mat(x).shuff5_rest(:,y);
    above95 = [above95, a > b]; %binary vector where CCG passed 95% confidence interval
    below5 = [below5, a < c]; %binary vector where CCG below 5% confidence interval
end

out(zz).delta = alignDelta_rest;
out(zz).delta50 = shuff50_rest;
out(zz).delta95 = shuff95_rest;
out(zz).align = align_rest;
out(zz).above95 = above95;
out(zz).below5 = below5;
out(zz).time = time;
end

end
