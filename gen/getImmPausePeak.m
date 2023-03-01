function [amp, dur, freq, thres] = getImmPausePeak(beh, varargin)
%Extract characteristics of troughs and peaks of ACh signal during periods
%of immobility
%
%   [amp, dur, freq, thres] = getImmPausePeak(beh)
%   [amp, dur, freq, thres] = getImmPausePeak(beh, NumStd)
%
%   Description: This function is analyzing fluctuations of ACh signal
%   during periods of immobility, and extracting the frequency, duration,
%   and amplitude of the troughs and peaks of fluctuations.
%       Frequency is defined as the number of events per second.
%       Duration is defined as the width at the half maximum amplitude
%       Amplitude is defined as the amplitude from zero
%
%   Input:
%   - beh - A data structure containing photometry, behavior
%   information for multiple recordings, using extractBeh or extractData
%   - optional inputs:
%       - NumStd - Number of standard deviations of filtered photometry
%       signal during immobility. Troughs and peaks are identified that
%       have an amplitude that is more than specified value * STD
%           default NumStd = 2
%
%   Output:
%   - dur - 2-column matrix with amplitudes for all troughs and peaks
%       amp(:,1) are amplitudes of TROUGHs
%       amp(:,2) are amplitudes of PEAKs
%   - dur - 2-column  matrix with duration (width at half max)
%       dur(:,1) are durations of TROUGHs
%       dur(:,2) for durations of PEAKs
%   - freq - 2-column  matrix with frequency
%       freq(:,1) are frequencies of TROUGHs
%       freq(:,2) for frequencies of PEAKs
%   - thres - Vector of thresholds (NumStd * STD) used for identification
%   of troughs and peaks
%
%   Author: Anya Krok, August 2022

    %% Initialize inputs
    thres = [];
    switch nargin
        case 1
            NumStd = 2; % default is 2
        case 2
            ii = varargin{1};
            if length(varargin{1}) > 1
                thres = ii; NumStd = 2;
            elseif length(varargin{1}) == 1
                NumStd = ii;
            end
    end
    
    %% Initialize outputs
    amp = nan(length(beh),2); % 1st column will contain amplitudes of troughs, 2nd column will contain amplitudes of peaks
    dur = nan(length(beh),2); 
    freq = nan(length(beh),2); 

    %% Analyze
    for x = 1:length(beh) % iterate over recordings
        if isempty(beh(x).Fs); continue; end
        Fs = beh(x).Fs;
        fp_mat = beh(x).FP{1}; 
        fp_mat = fp_mat - nanmean(fp_mat);

        rewWindow = Fs; % how many samples after reward delivery is the reward window
        idx_rew = [];
        if isfield(beh, 'reward'); if ~isempty(beh(x).reward)
            idx_rew = extractEventST([1:length(fp_mat)]', floor(beh(x).reward), floor(beh(x).reward)+rewWindow, 1); % identify sample during reward
            end; end
        idx_rest = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % index: infusion + rest
        if isempty(idx_rest); continue; end
        idx_rest = idx_rest(~ismember(idx_rest, idx_rew)); % exclude reward, include rest

        %% Bandpass filter
        Fpass = [0.5 4];
        % Fs = 50; %sampling rate, has to be at least double of your high pass frequency
        Wn = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; 
        [b,a] = butter(3,Wn);
        data_filt = filtfilt(b, a, fp_mat(:,1)); % signal is your photometry data, output is the filtered data

        %% Instantaneous phase
        data_filt = data_filt(idx_rest); % IMMOBILITY ONLY
        H = hilbert(double(data_filt));
        data_phase  = angle(H); % output is the instantaneous phase
        fp_phase = data_phase;
        fp_deg = rad2deg(data_phase);

        %Take the absolute of the filtered signal and calculate the standard deviation
        rmssig  = abs(data_filt);
        stdsig = std(rmssig);        

        %Find the index of peaks that are more than NumStd standard deviations (NumStd could be 1.5 standard deviations)
        % finds the 0 degree phase indices that cross 1.5 standard deviations. Might have to remove the ; depending on whether your data is a row or column array
        if length(varargin{1}) == 1
            thres(x) = NumStd*stdsig;
        end
        idxPeak = find(data_filt>thres(x) & [0; diff(data_phase>0)]); 
        idxPause = find(-data_filt>thres(x) & [0; diff(-data_phase>0)]); 
    %         figure; hold on
    %         plot(data_filt, 'k')
    %         stem(idxPeak, 10*ones(length(idxPeak),1), 'g'); 
    %         stem(idxPause, 10*ones(length(idxPause),1), 'r');

        %% Peak characterization  
        tmp_amp_peak = data_filt(idxPeak); 
        tmp_amp_pause = data_filt(idxPause);
        tmp_dur_peak = []; tmp_dur_pause = []; win = 1*Fs;
        for z = 1:length(idxPeak)
            halfMax = 0.5*data_filt(idxPeak(z)); % amplitude at half-max
            if idxPeak(z)-win <= 0
                a = [data_filt(1 : idxPeak(z))]; a = flipud(a); % if peak is close to start of recording
            else
                a = [data_filt(idxPeak(z)-win : idxPeak(z))]; a = flipud(a); % segment preceding idx of maximum deflection
            end
            j2 = find(a < 0, 1, 'first') - 1; 
            if ~isempty(j2); a = a(1:j2); end % retain only values above 0
            c = []; [~,c(1)] = min(abs(a-halfMax)); c(1) = c(1)-1; % idx at half max
            if idxPeak(z)+win > length(data_filt)
                a = [data_filt(idxPeak(z) : end)]; 
            else
                a = [data_filt(idxPeak(z) : idxPeak(z)+win)]; % segment following idx of maximum deflection
            end

            j2 = find(a < 0, 1, 'first') - 1;
            if ~isempty(j2); a = a(1:j2); end % retain only values above 0
            [~,c(2)] = min(abs(a-halfMax)); % idx at half max

            tmp_dur_peak(z) = sum(c);
        end % width at half max for peaks
        for z = 1:length(idxPause)
            halfMax = 0.5*data_filt(idxPause(z)); % amplitude at half-max
            if idxPause(z)-win <= 0
                a = [data_filt(1 : idxPause(z))]; a = flipud(a); % if peak is close to start of recording
            else
                a = [data_filt(idxPause(z)-win : idxPause(z))]; a = flipud(a); % segment preceding idx of maximum deflection
            end
            j2 = find(a > 0, 1, 'first') - 1;
            if ~isempty(j2); a = a(1:j2); end % retain only values above 0
            c = []; [~,c(1)] = min(abs(a-halfMax)); c(1) = c(1)-1; % idx at half max
            if idxPause(z)+win > length(data_filt)
                a = [data_filt(idxPause(z) : end)]; 
            else
                a = [data_filt(idxPause(z) : idxPause(z)+win)]; % segment following idx of maximum deflection
            end
            j2 = find(a > 0, 1, 'first') - 1;
            if ~isempty(j2); a = a(1:j2); end % retain only values above 0
            [~,c(2)] = min(abs(a-halfMax)); % idx at half max

            tmp_dur_pause(z) = sum(c);
        end % width at half max for pauses

        %%
        amp(x,1) = nanmean(tmp_amp_pause); % Amplitude of maximum deflection
        amp(x,2) = nanmean(tmp_amp_peak);
        dur(x,1) = nanmean(tmp_dur_pause); % Duration at half maximum
        dur(x,2) = nanmean(tmp_dur_peak); 
        freq(x,1) = (1./nanmean(diff(idxPause)))*Fs; % Frequency of maximum deflections
        freq(x,2) = (1./nanmean(diff(idxPeak)))*Fs; % Frequency of maximum deflections
    end
    
end