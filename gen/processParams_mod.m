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
params.acqFs = 5000; % Acquisition sampling frequency
params.dsType = 2; % 1 = Bin Summing; 2 = Bin Averaging; 3 = Traditional (NOT RECOMMENDED)
params.dsRate = 100; % Downsampling rate if you want to downsample the signal
%This dsRate will also be applied to all signals during the analysis
%pipeline
%dsRate = 1 if you do not want to downsample

%% Filter Parameters
params.FP.lpCut = 10; % Cut-off frequency for filter
params.FP.filtOrder = 8; % Order of the filter

%% Baseline Parameters
params.FP.basePrc = 5; % Percentile value from 1 - 100 to use when finding baseline points
%Note: Lower percentiles are used because the mean of signal is not true
%baseline
params.FP.winSize = 10; % Window size for baselining in seconds
params.FP.winOv = 0; %Window overlap size in seconds
params.FP.interpType = 'linear'; % 'linear' 'spline' 
params.FP.fitType = 'interp'; % Fit method 'interp' , 'exp' , 'line'

%% Demodulation Parameters
%When demodulating signals, the filter creates edge artifacts. We record
%for a few seconds longer, so we can remove x seconds from the beginning
%and end
%Adjust the variable to "0" if it's a normal photometry recording
params.FP.sigEdge = 15; %Time in seconds of data to be removed from beginning and end of signal

params.FP.modFreq = [217 319];

%% Behavior Parameters Parameters
%Wheel Parameters
params.mov.radius = 9.8; %Radius of the wheel used. Note it can be meters or centimeters. Just keep track of your units
params.mov.winSize = 0.5; %This is the window size for the moving avg filter applied to unwrapped encoder data 500ms windows work well
%Onset/Offset Parameters
%Movement Onset and Offset Parameters
params.mov.velThres = 4; %(same units as radius)/s
params.mov.minRunTime = 4; %Threshold for minimum time spent running for movement bouts (in seconds)
params.mov.minRestTime = 4; %Threshold for minimum time spent rest for movement bout (in seconds)
params.mov.finalOnset = 0; %Boolean value -- Decides if you want to include or exlcude the final 
% onset if the acquisition ends before the offset
params.mov.timeThres = 4; %Make sure a bout is above a certain time-length
params.mov.timeBefore = 4; %Time to display preceding movement onset and offset
params.mov.timeAfter = 4; %Time to display following movement onset and offset
params.mov.iterSTD = 0.5; %Minimum iteration std value
params.mov.iterWin = 3; %Window size used to find minimum iteration value
%Rest Onset and Offset Parameters
params.mov.minRestTime_rest = 4;
params.mov.minRunTime_rest = 1;
params.mov.velThres_rest = 0.25;
params.mov.timeThres_rest = 4;
params.mov.timeShift_rest = 0.5;

%% Peak Analysis
% params.peaks.minHeight = 0.4;
% params.peaks.minProminence = 0.4;
% params.peaks.smoothWin = 2;
% params.troughs.minHeight = 0.4;
% params.troughs.minProminence = 0.4;

%% Cross-Correlations
%params.cc.lag = 1; %In seconds how much to shift forward and backwards
%params.cc.type = 1; % 1 = Cross-Correlation 2 - Cross-Covariance

%% Opto-Pulse Analysis
%   'threshold' - a.u. or V, depends on voltage output of pulse generator
%       %arduino(for in vivo): 4V
%       %wavesurfer(photometry): 0.15V
params.opto.threshold = 1; 
params.opto.stimtype = 'excitation';
params.opto.cutoff = 20; %cutoff freq 20Hz
params.opto.order = 10;
params.opto.filtType = 'lowpass';
params.opto.dsRate = params.dsRate; 
params.opto.dsType = params.dsType;

%%
answer = inputdlg({...
    'Sampling freq (Hz)','Downsampling Rate',...
    'Behavior type (open/wheel)',...
    'Velocity threshold for wheel','Rest threshold for wheel'...
    'Photometry? (yes/no)','Optogenetics? (yes/no)','Wheel Radius (cm)'},...
    'Input', 1,...
    {num2str(params.acqFs),num2str(params.dsRate),'wheel',...
    num2str(params.mov.velThres),num2str(params.mov.velThres_rest),...
    'yes','no',num2str(params.mov.radius)});

    params.dsRate = str2num(answer{2}); %Downsampling Rate
    
    if ~strcmp(answer{3},'open') && ~strcmp(answer{3},'wheel')
        error('params_mov did not have proper input for behavior type');
        return
    end
    switch answer{3} %Behavior type (corridor/wheel)
        case 'open'
            params.FP.fitType = 'interp'; %interp if freely moving
            params = rmfield(params,'mov'); 
        case 'wheel'
            params.FP.fitType = 'interp'; %interp if head-fixed
            params.mov.velThres = str2num(answer{4});
            params.mov.velThres_rest = str2num(answer{5});
            params.mov.radius = str2num(answer{8});
    end
    
    if ~strcmp(answer{6},'yes') && ~strcmp(answer{6},'no')
        error('params_mov did not have proper input for photometry yes/no');
        return
    end
    switch answer{6} %Photometry? (yes/no)
        case 'yes'
            answerFP = inputdlg({'FP fitType: interp',...
                'FP modulation? (yes/no)',...
                'FP sigEdge: 0 for none',...
                'modFreq-1 (Hz)','modFreq-2 (or isoFreq) (Hz)'},'FP',1,...
                {params.FP.fitType,'yes',num2str(params.FP.sigEdge),...
                num2str(params.FP.modFreq(1)),num2str(params.FP.modFreq(2))});
            params.FP.fitType = answerFP{1};
            switch answerFP{2}
                case 'yes'
                    params.FP.sigEdge = str2num(answerFP{3});
                    params.FP.modFreq = [str2num(answerFP{4}), str2num(answerFP{5})];
                case 'no'
                    params.FP.sigEdge = 0;
                    params.FP = rmfield(params.FP,{'modFreq'});
            end
        case 'no'
            params = rmfield(params,'FP'); params.FP.sigEdge = 0;
            params = rmfield(params,'peaks');
    end
    
    if ~strcmp(answer{7},'yes') && ~strcmp(answer{7},'no')
        error('params_mov did not have proper input for optogenetics yes/no');
        return
    end
    switch answer{7} %Optogenetics? (yes/no)
        case 'yes'
            answerOp = inputdlg({'Opto stimtype:','Opto dsRate:','Stim freq (Hz):',...
                '#pulse/train:','Stim ISI:'},...
                'Opto',1,...
                {params.opto.stimtype,num2str(params.dsRate),'','',''});
            params.opto.stimtype = answerOp{1};
            params.opto.dsRate = str2num(answerOp{2});
            params.opto.freq = str2num(answerOp{3});
            params.opto.train = str2num(answerOp{4});
            params.opto.ISI = answerOp{5};
        case 'no'
            params = rmfield(params,'opto');
    end
