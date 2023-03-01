function [spikeTimes, spikeMat, tVec] = poissonSpikeGen(fr, tSim, nTrials)
% Poisson process to simulate spike train
%
% [spikeTimes, spikeMat, tVec] = poissonSpikeGen(fr, tSim, nTrials)
%
% Description: This function uses Poisson process to 
% simulate spike trains that have characteristics close to real neurons. 
%
% INPUTS
%   fr - firing rate
%   tSim - duration of simulation
%   nTrials - number of trials
%
% OUTPUTS
%   spikeTimes - cell array of spike times, in seconds
%   spikeMat - binary spike matrix
%   tVec - time vector
%       In time vector tVec, each time stamp is 1 ms long (dt). 
%       The last time stamp will start at tSim - dt, which represents the time interval between tSim - dt and tSim.
%
% Anya Krok, March 2021

dt = 1/1000; % Set dt = 1 ms
nBins = floor(tSim/dt);
spikeMat = rand(nTrials , nBins) < fr*dt;
tVec = 0:dt:tSim-dt;
spikeTimes = {};
for x = 1:nTrials
    spikeTimes{x} = find(spikeMat(x,:)).*dt; 
end

% 1. Compute the product of fr*dt. 
%   Let?s stimulate a 10 ms long spike train for a neuron firing at 100 Hz. 
%   Thus, fr = 100 Hz, dt = 1 ms and fr*dt = 0.1 (remember that Hz is the 
%   inverse of s and 1 ms is 1/1000 s, so fr*dt is dimensionless).
%
% 2. Generate uniformly distributed random numbers between 0 and 1. 
%   In MATLAB, use the function rand() to do so. rand(1, 10) will generate 10 random numbers.
%
% 3. Compare each random number to fr*dt. 
%   If the product is less than fr*dt (x < fr*dt), then there is a spike.
end
