function beh = extractBeh(varargin)
%%Extract processed data from multiple recordings into single structure
%
% [beh] = extractBeh() - load from files, will be prompted to select
% [beh] = extractBeh(fPath, fName) - load from specified files
% [beh] = extractBeh(beh) - to add to existing structure
% [beh] = extractBeh(fPath, fName, beh) - to add to existing structure
%
% Description: Extract processed photometry, locomotion, and/or reward data 
% from multiple individual recording files into a larger structure to be 
% used for data analysis and plotting
%
% INPUTS
%   'fPath' - Character array containing folder path where data files are
%       example: 'R:\tritsn01labspace\Anya\FiberPhotometry\AK201-206\220105'
%   'fName' - Cell array, with each cell containing file names for each
%   recording to be added to structure
%   'beh' (optional) - Structure previously created with extractBeh
%
% OUPUTS
%   'beh' - Structure with data from multiple recordings
%
% Anya Krok, January 2021

    %% INPUTS
    switch nargin
        case 0
            [fName,fPath] = uigetfile('*.mat','Select data files to add to beh structure','MultiSelect','On');
            beh = struct;
        case 1
            [fName,fPath] = uigetfile('*.mat','Select data files to add to beh structure','MultiSelect','On');
            beh = varargin{1};
        case 2
            fPath = varargin{1};
            fName = varargin{2};
            beh   = struct;
        case 3
            fPath = varargin{1};
            fName = varargin{2};
            beh   = varargin{3};
    end
    if ~iscell(fName); fName = {fName}; end
    fName = sort(fName);
    
    %%
    for f = 1:length(fName) 
        load(fullfile(fPath,fName{f})); 
        [an,b] = strtok(fName{f},'_'); day = strtok(b,'_');
        x = 1+length(beh);
        beh(x).rec = [an,'-',day]; 
        beh(x).site = 'DLS'; % CHANGE, or adjust later
        beh(x).task = 'wheel';

        %% PHOTOMETRY
        beh(x).Fs = data.gen.Fs; % Sampling frequency, in Hz
        beh(x).time = data.final.time; % Time vector
        if isfield(data.final,'FP')
            beh(x).FP = data.final.FP;
            beh(x).nbFP = data.final.nbFP; % Photometry signal(s)
            beh(x).FPnames = data.final.FPnames; % Names of photometry signal(s)
        end

        %% LOCOMOTION
        if isfield(data.final,'vel') % If movement data exists
            [~,ii] = max(abs(data.final.vel));
            if data.final.vel(ii) < 0
                data.final.vel = -data.final.vel; % Flip velocity for recordings on IV rig#1, which has inverted positional encoder signal
                fprintf('%s - inverted velocity signal. Flipping and re-saving \n',beh(x).rec);
                save(fullfile(fPath,fName{f}),'data'); % Overwrite data file to adjust
            end
            beh(x).vel = data.final.vel; % Velocity signal
            if isfield(data.final,'mov')
                if isfield(data.final.mov,'onsets')
                    beh(x).on = data.final.mov.onsets; 
                    beh(x).off = data.final.mov.offsets; % Locomotion onset/offset times in sampling freq (data.gen.Fs), NOT in seconds
                else
                    beh(x).on = []; beh(x).off = [];
                end
                if isfield(data.final.mov,'onsetsRest')
                    beh(x).onRest = data.final.mov.onsetsRest; 
                    beh(x).offRest = data.final.mov.offsetsRest; % Immobility onset/offset times in sampling freq (data.gen.Fs), NOT in seconds
                else
                    beh(x).onRest = []; beh(x).offRest = [];
                end
            end
        end

        %% REWARD
        if isfield(data.acq,'rew')
            beh(x).task = 'reward';
            beh(x).reward = data.final.rew.onset;    % Reward delivery time in sampling freq (data.gen.Fs), NOT in seconds
            beh(x).lick = data.final.lick.onset;     % Lick times in sampling freq (data.gen.Fs), NOT in seconds
        end
        
        %% CAMERA
        if isfield(data.acq,'cam')
            beh(x).task = 'openField';
            beh(x).cam = data.final.cam.on;     % Camera trigger times in sampling freq (data.gen.Fs), NOT in seconds
        end

        %%
        fprintf('Extracted from %s\n',fName{f});
    end
    if isempty(beh(1).Fs); beh(1) = []; end
end