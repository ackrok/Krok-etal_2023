function stSub = extractEventST(st,on,off,flag)
% Extract spike times within a specified event time range
%
% stSub = extractEventST(st,on,off)
%
% Description: This function will take inputs of spike times and the start
% and end times for an event, and return a vector of spike times that occur
% within the event time range.
%
% INPUT
%   'st' - vector of spike times, in units that match that of event times
%   'on' - vector of start times for event(s)
%   'off' - vector of end times for event(s)
%   'flag' - 0 or 1.
%       0: concatenate ISIs, then regenerate STs
%       1: concatentae STs
%
% OUTPUT
%   'stSub' - vector of subset of spike times occuring within event range
%
% Anya Krok, April 2020
% Updated: April 26 2020, adjust spike times for start time of each event,
% in order to limit large ISIs between spikes of adjacent events
% Updated: March 3 2021, added flag

if length(on) ~= length(off) 
    error('Start and end times for event(s) do not match');
end

switch flag
    case 0
        ISI = []; 
        for x = 1:length(on)
            range = [on(x), off(x)];                 % Sample range
            check = st >= range(1) & st <= range(2); % Check if values within range
            stTmp = st;
            stTmp(~check) = nan; % If no values within range, return NaN
            stTmp = stTmp(~isnan(stTmp)); % Remove NaNs from spike time vector
            ISI = [ISI; diff(stTmp)]; % concatenate ISI between spike times within event range
        end
        stSub = cumsum(ISI); % Cumulative sum of ISIs returns spike time 

    case 1
        stSub = []; 
        for x = 1:length(on)
            range = [on(x), off(x)];                 % Sample range
            check = st >= range(1) & st <= range(2); % Check if values within range
            stTmp = st;
            stTmp(~check) = nan; % If no values within range, return NaN
            stTmp = stTmp(~isnan(stTmp)); % Remove NaNs from spike time vector
            stSub = [stSub; stTmp]; % concatenate ISI between spike times within event range
        end

end
end