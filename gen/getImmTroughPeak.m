function cannula = getImmTroughPeak(cannula)
%% 
% cannula = getImmTroughPeak(cannula)
%
% Anya Krok, March 2022

%% INPUTS
NumStd = 1.5;

%%
params = [];
for y = 1:length(cannula)
    %%
    beh = cannula(y).s; % extract beh sub-structure within cannula structure
    amp = nan(length(beh),2); % initialize output matrices
    dur = nan(length(beh),2);  % initialize output matrices
    freq = nan(length(beh),2); % initialize output matrices
    ach2ach = cell(length(beh),2); % initialize output matrices
    da2ach = cell(length(beh),2); % initialize output matrices

    idxStates = extractBehavioralStates(beh); % extract indices for each behavioral state
    for x = 1:length(beh) % iterate over animal
        if isempty(beh(x).Fs); continue; end
        fpMat = [beh(x).FP{1}, beh(x).FP{2}];
        fpMat = fpMat - nanmean(fpMat); % subtract mean of recording to center around 0%
        Fs = beh(x).Fs;
        idxImm = idxStates{x,1}; % extract samples during immobility, excluding reward and locomotion
        winInf = [cannula(y).win(x,1).*(Fs*60), cannula(y).win(x,2).*(Fs*60)]'; % infusion window
        idxImmInf = idxImm(idxImm > winInf(1) & idxImm < winInf(2)); % immobility during infusion window

        %% Bandpass filter
        Fpass = [0.1 10];
        Fs = beh(x).Fs; %sampling rate, has to be at least double of your high pass frequency
        Wn = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; 
        [b,a] = butter(3,Wn);
        dataFilt = filtfilt(b,a,fpMat(:,1)); % signal is your photometry data, output is the filtered data

        %% Instantaneous phase
        dataFilt = dataFilt(idxImmInf); % IMMOBILITY ONLY
        H = hilbert(double(dataFilt));
        dataPhase  = angle(H); % output is the instantaneous phase
%         fpPhase = dataPhase;
%         fpDeg = rad2deg(dataPhase);

        %Take the absolute of the filtered signal and calculate the standard deviation
        rmssig  = abs(dataFilt);
        stdsig = std(rmssig);        

        %Find the index of peaks that are more than NumStd standard deviations (NumStd could be 1.5 standard deviations)
        % finds the 0 degree phase indices that cross 1.5 standard deviations. Might have to remove the ; depending on whether your data is a row or column array
%         NumStd = 1.5;
        % if y == 1; params(x) = NumStd*stdsig; end % Use threshold for saline infusion, which should be in cannula(1)
        params(x) = NumStd*stdsig;
        idxPeak = find(dataFilt>(params(x)) & [0; diff(dataPhase>0)]); 
        idxTrgh = find(-dataFilt>(params(x)) & [0; diff(-dataPhase>0)]); 
%         figure; hold on
%         plot(dataFilt, 'k')
%         stem(idxPeak, 10*ones(length(idxPeak),1), 'g'); 
%         stem(idxTrgh, 10*ones(length(idxPause),1), 'r');

        %% Peak characterization
        tmpDurPeak = []; tmpDurTrgh = []; win = 1*Fs;
        for z = 1:length(idxPeak)
            halfMax = 0.5*dataFilt(idxPeak(z)); % amplitude at half-max
            if idxPeak(z)-win < 0
                a = [dataFilt(1 : idxPeak(z))]; a = flipud(a); % if peak is close to start of recording
            else
                a = [dataFilt(idxPeak(z)-win : idxPeak(z))]; a = flipud(a); % segment preceding idx of maximum deflection
            end
            b = find(a < 0, 1, 'first') - 1; 
            if ~isempty(b); a = a(1:b); end % retain only values above 0
            c = []; [~,c(1)] = min(abs(a-halfMax)); c(1) = c(1)-1; % idx at half max
            if idxPeak(z)+win >= length(dataFilt)
                a = [dataFilt(idxPeak(z) : end)]; 
            else
                a = [dataFilt(idxPeak(z) : idxPeak(z)+win)]; % segment following idx of maximum deflection
            end
                
            b = find(a < 0, 1, 'first') - 1;
            if ~isempty(b); a = a(1:b); end % retain only values above 0
            [~,c(2)] = min(abs(a-halfMax)); % idx at half max
            tmpDurPeak(z) = sum(c);
        end % width at half max for peaks
        for z = 1:length(idxTrgh)
            halfMax = 0.5*dataFilt(idxTrgh(z)); % amplitude at half-max
            if idxTrgh(z)-win <= 0
                a = [dataFilt(1 : idxTrgh(z))]; a = flipud(a); % if peak is close to start of recording
            else
                a = [dataFilt(idxTrgh(z)-win : idxTrgh(z))]; a = flipud(a); % segment preceding idx of maximum deflection
            end
            b = find(a > 0, 1, 'first') - 1;
            if ~isempty(b); a = a(1:b); end % retain only values above 0
            c = []; [~,c(1)] = min(abs(a-halfMax)); c(1) = c(1)-1; % idx at half max
            if idxTrgh(z)+win > length(dataFilt)
                a = [dataFilt(idxTrgh(z) : end)]; 
            else
                a = [dataFilt(idxTrgh(z) : idxTrgh(z)+win)]; % segment following idx of maximum deflection
            end
            b = find(a > 0, 1, 'first') - 1;
            if ~isempty(b); a = a(1:b); end % retain only values above 0
            [~,c(2)] = min(abs(a-halfMax)); % idx at half max
            tmpDurTrgh(z) = sum(c);
        end % width at half max for troughs

        amp(x,1) = nanmean(dataFilt(idxTrgh)); % Amplitude of maximum deflection for troughs
        amp(x,2) = nanmean(dataFilt(idxPeak)); % Amplitude of maximum deflection for peaks
        dur(x,1) = nanmean(tmpDurTrgh).*(1000/Fs); % Duration at half max for troughs, adjust from samples to ms
        dur(x,2) = nanmean(tmpDurPeak).*(1000/Fs); % Duration at half max for peaks, adjust from samples to ms
        freq(x,1) = (1./nanmean(diff(idxTrgh)))*Fs; % Frequency of maximum deflections for troughs
        freq(x,2) = (1./nanmean(diff(idxPeak)))*Fs; % Frequency of maximum deflections for peaks
        
        %% Photometry to peaks/pauses
        [sta, t] = getSTA(fpMat(idxImmInf,1), idxTrgh/Fs, Fs, [-6 2]); % ACh photometry (or whichever signal is in beh(x).FP{1}) aligned to troughs for ACh signal
        ach2ach{x,1} = sta;
        sta = getSTA(fpMat(idxImmInf,1), idxPeak/Fs, Fs, [-6 2]); % ACh photometry (or whichever signal is in beh(x).FP{1}) aligned to peaks for ACh signal
        ach2ach{x,2} = sta;
        sta = getSTA(fpMat(idxImmInf,2), idxTrgh/Fs, Fs, [-6 2]); % DA photometry (or whichever signal is in beh(x).FP{2}) aligned to troughs for ACh signal
        da2ach{x,1} = sta;
        sta = getSTA(fpMat(idxImmInf,2), idxPeak/Fs, Fs, [-6 2]); % DA photometry (or whichever signal is in beh(x).FP{2}) aligned to peaks for ACh signal
        da2ach{x,2} = sta;
    end
    cannula(y).params = params;
    cannula(y).amp = amp; 
    cannula(y).dur = dur; 
    cannula(y).freq = freq;
    cannula(y).ach2ach = ach2ach;
    cannula(y).da2ach = da2ach;
    cannula(y).staTime = t;
end
fprintf('Analysis of ACh troughs and peaks during immobility done.\n');
end