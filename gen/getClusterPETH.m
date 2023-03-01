function [peth] = getClusterPETH (st, events, bin, window)
% Peri-Event Time Histogram (PETH) - alignment of spike times to event time 
%
% [peth] = getClusterPETH (st, events, bin, window)
%
% INPUT
% Note: all inputs must be in same time scale (e.g. seconds)
%   'st' - cell array with spike times
%   'events'  -  vector with event times
%   'bin' - bin width (s) for binning of spike times
%   'window' - window (s) to analyze around event time stamp
%
% OUTPUT
%   'peth' - structure with PETH output
%       'peth.fr' - peri-event time histogram bin counts
%                   adjusted for bin width and number of events
%       'peth.fr_z' - rudimentary z-scored peri event rate
%       'peth.time' - time vector for plotting
%
% Created by: Anya Krok, February 2019
% Updated: Anya Krok, March 2020
%

%% Input Variables
% bin = peth.bin; %Bin width for binning of spike times
% window = peth.window; %Window to analyze around event time stamp
edges = [window(1): bin :window(2)]; edges = edges(:); %Edges of bins within which spike times were binned for PETH
if ~iscell(st); st = {st}; end

%% Compute PETH
cts_byUnit = zeros(length(edges)-1,length(st)); %Initialize temporary matrix for storing average counts for each unit
cts_byUnit(cts_byUnit == 0) = NaN; %Convert to NaN matrix
save_cts_byEvent = cell(length(st),1); 
for x = 1:length(st) %Iterate over units
    cts_byEvent = zeros(length(edges)-1,length(events)); %Initialize temporary matrix for storing counts for each event
    cts_byEvent(cts_byEvent == 0) = NaN;
    if ~isempty(st{x})
        for y = 1:length(events) %Iterate over events
            nCounts = histcounts(st{x}, ... %Compute PETH - bin spike times within range surrounding event time
                'BinWidth', bin, ...                %Specify width of bins 
                'BinLimits', events(y) + window);   %Specify range surrounding event time
            cts_byEvent (:,y) = nCounts./bin;     %Adjust for width of each bin
        end
        cts_byUnit (:,x) = nansum(cts_byEvent,2)./length(events); %Adjust for number of events
    end
    save_cts_byEvent{x} = cts_byEvent.*bin;
end

%% Compute z-score for firing rate
%[fr_z, mu, sigma] = zscore(cts_byUnit,[],1);
mu = []; sigma = []; fr_z = [];
for x = 1:length(st)
    mu = 1/mean(diff(st{x})); sigma = std(diff(st{x}));
    fr_z(:,x) = (cts_byUnit(:,x) - mu)./sigma;
end

%% Save Outputs
time = edges + 0.5*bin; %To convert edges into time vector, shift half bin forward
time = time(1:end-1); %and remove last value - resulting vector t = 0 falls in between bins instead of at center of bin

peth.time = time; %Time vector for plotting
peth.fr   = cts_byUnit; %PETH by unit
peth.fr_z = fr_z; %PETH z-scored for each unit
peth.cts  = save_cts_byEvent; %Counts, adjusted by bin size

end
