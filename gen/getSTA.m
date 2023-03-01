function [sta, time, sta_z] = getSTA(signal, events, Fs, window)
% Spike Triggered Average (STA) - alignment of signal to event time 
%
% [sta, time, sta_z] = getSTA(signal, events, Fs, window)
%
% example: getSTA(fp, beh_onsets, 50, [0 5]); to extract photometry signal
% from t = 0s after movement onset to t = 5s after movement onset
%
% INPUT
%   'signal' - vector with signal that you wish to align to events, in
%   sampling frequency defined by Fs
%   'Fs' - sampling frequency of signal
%   'events'  -  vector with event times, in seconds
%   'window' - window to analyze around events, in seconds
%
% OUTPUT
%   'sta' - matrix with event-aligned signal
%       matrix will have numRows = Fs*[window(2)-window(1) 
%       and numColumns = numEvents
%   'time' - time vector for plotting
%   'sta_z' - matrix with rudimentary z-scored event-aligned signal
%
% Created by: Anya Krok, February 2019
% Updated: Anya Krok, March 2020
%

range = [window(1)*Fs : window(2)*Fs]; nRange = length(range); % Window range, in samples, matching signal
time  = range(:)./Fs; % Time vector, in seconds

sta = zeros(nRange, length(events)); % Initialize matrix with out
mu = mean(signal); sigma = std(signal); % Compute mu, sigma over entire signal

for x = 1:length(events)
    eventFs = round(events(x).*Fs); % Event time matching sampling rate of singal
    if eventFs + range(1) >= 1 && eventFs + range(end) <= length(signal)
        sta(:,x) = signal(eventFs + range); % Extract signal from window around this event
    else
        sta(:,x) = nan(nRange, 1); % If window extends beyond range of signal, fill with NaNs
    end
    sta_z(:,x) = (sta(:,x) - mu)./sigma; % z-score STA using mu and std over entire signal
end
