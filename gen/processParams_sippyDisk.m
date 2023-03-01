%% Process Parameters for Photometry Analysis
%Created By: Pratik Mistry
%Created On: 31 January 2019
%Edited On: 26 September 2019
%
%Description: This is a script with different variables whose values can be
%adjusted depending on the photometry signal that is being processed. The
%name of this file can and should be changed depending on the method and
%GECI used. Please read comments associated with variable
%
%
%% General Parameters
params.acqFs = 20000; % Acquisition sampling frequency
params.dsType = 2; % 1 = Bin Summing; 2 = Bin Averaging; 3 = Traditional (NOT RECOMMENDED)
params.dsRate = 400; % Downsampling rate if you want to downsample the signal
%This dsRate will also be applied to all signals during the analysis
%pipeline
%dsRate = 1 if you do not want to downsample

%% Filter Parameters
params.FP.lpCut = 15; % Cut-off frequency for filter
params.FP.filtOrder = 8; % Order of the filter

%% Baseline Parameters
params.FP.basePrc = 5; % Percentile value from 1 - 100 to use when finding baseline points
%Note: Lower percentiles are used because the mean of signal is not true
%baseline
params.FP.winSize = 5; % Window size for baselining in seconds
params.FP.winOv = 0; %Window overlap size in seconds
params.FP.interpType = 'linear'; % 'linear' 'spline' 
params.FP.fitType = 'interp'; % Fit method 'interp' , 'exp' , 'line'

%% Demodulation Parameters
%When demodulating signals, the filter creates edge artifacts. We record
%for a few seconds longer, so we can remove x seconds from the beginning
%and end
%Adjust the variable to "0" if it's a normal photometry recording
params.FP.sigEdge = 0; %Time in seconds of data to be removed from beginning and end of signal
% params.FP.modFreq = [217 319];

%% Behavior Parameters Parameters
%Wheel Parameters
params.mov.radius = 9.8; %Radius of the wheel used. Note it can be meters or centimeters. Just keep track of your units
params.mov.winSize = 0.25; %This is the window size for the moving avg filter applied to unwrapped encoder data 500ms windows work well
%Onset/Offset Parameters
%Movement Onset and Offset Parameters
params.mov.velThres = 0.001; %(same units as radius)/s
params.mov.minRunTime = 1; %Threshold for minimum time spent running for movement bouts (in seconds)
params.mov.minRestTime = 4; %Threshold for minimum time spent rest for movement bout (in seconds)
params.mov.finalOnset = 0; %Boolean value -- Decides if you want to include or exlcude the final 
% onset if the acquisition ends before the offset
params.mov.timeThres = 0.25; %Make sure a bout is above a certain time-length
params.mov.timeBefore = 2; %Time to display preceding movement onset and offset
params.mov.timeAfter = 2; %Time to display following movement onset and offset
params.mov.iterSTD = 0.5; %Minimum iteration std value
params.mov.iterWin = 1; %Window size used to find minimum iteration value
%Rest Onset and Offset Parameters
params.mov.minRestTime_rest = 1;
params.mov.minRunTime_rest = 1;
params.mov.velThres_rest = 0.5;
params.mov.timeThres_rest = 2;
params.mov.timeShift_rest = 0.5;
