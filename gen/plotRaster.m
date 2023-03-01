function [fig] = plotRaster(spikeMat, tVec)
% Plot the spike matrix generated using poissonSpikeGen
%
% [] = plotRaster(spikeMat, tVec)
%
% INPUTS
%   spikeMat - spike matrix, from poissonSpikeGen
%   tVec - time vector, from poissonSpikeGen
%
% OUTPUTS
%   fig - figure handle
%
% Anya Krok, March 2021

fig = figure;
hold all;
for trialCount = 1:size(spikeMat,1)
    spikePos = tVec(spikeMat(trialCount, :));
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4], 'k');
    end
end
ylim([0 size(spikeMat, 1)+1]);
