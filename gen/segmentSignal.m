function [seg_fp, seg_samp] = segmentSignal(fp,winSize,winOv,Fs,varargin)
%segmentSignal - Extract segments of specified windows from photometry
%signal
%
%   [seg_fp] = segmentSignal(fp,winSize,winOv,Fs)
%   [seg_fp, seg_samp] = segmentSignal(fp,winSize,winOv,Fs,evnts)
%
%   Description: This code will use a moving window to extract segments of 
%   a continuous signal that are a specified length of time
%
%   Input:
%   - FP - Photometry signal to extract segments from
%   - winSize - Window size in seconds for finding baseline
%   - winOv - Overlap size in seconds for finding baseline
%   - Fs - Sampling Rate
%   - events (optional) - Two column matrix of event start and end times
%       , in samples
%
%   Output:
%   - seg_fp - Matrix of photometry signal segments
%   - seg_samp - Two columns of start and end samples of each segment
%
%   Author: Anya Krok, December 2021
%

%% Inputs
if size(fp,1) == 1
    fp = fp'; %Ensure the photomtery vector is a column vector
end

if nargin > 4
    events = varargin{1}; %Two columns of event start and end times, in samples
else
    events = [1, length(fp)];
end

%% Window
winSize = winSize * Fs; %Convert window size from seconds to samples
winOv = winOv * Fs;     %Convert overlap window from seconds to samples
%Determine the step size of the window:
%If the overlap is 0 or empty then it will use the window size as the step
%size. If the overlap is greater than 0 the step size will be the window
%size subtracted by the overlap size
if winOv == 0 || isempty(winOv)
    winStep = winSize;
else
    winStep = winSize - winOv;
end

%% Segment photometry
%fp_seg is a matrix of zeros that will contain all concatenated segments
seg_fp = []; seg_samp = [];

for e = 1:size(events,1) %Iterate over number of events
    fp_ev = fp(events(e,1):events(e,2));
    Ls = length(fp_ev);        %Get length of photometry trace
    L = 1:Ls; L = L';       %Create a column vector from 1 to total data points in trace
    nSeg = floor(Ls/(winSize-winOv)); %Determine number of segments

    %X is a matrix of start and end samples of segments
    %Y is a matrix of zeros that will contain segments
    X = zeros(nSeg,2); 
    Y = zeros(winSize+1,nSeg);

    %The following for loop goes through the photometry vector and finds
    %baseline values of the windowed photometry trace according to a certain
    %percentile
    for n = 0:nSeg-1
        I1 = (n*winStep)+1;
        I2 = I1 + winSize;
        if I2>Ls
            I2 = Ls; tmp = fp_ev(I1:I2); 
            tmp = [tmp;nan(size(Y,1)-length(tmp),1)]; % Append nans to end of segment when shorter than window
            Y(:,n+1) = tmp;
        else
            Y(:,n+1) = fp_ev(I1:I2); % Store photometry segment into matrix
        end
        X(n+1,:) = [I1 I2] + events(e,1) - 1;
    end

    %Append to fp_seg matrix
    seg_fp = [seg_fp, Y]; seg_samp = [seg_samp; X];
end
end