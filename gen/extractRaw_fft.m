function [rawS] = extractRaw_fft(varargin)
%%Extract raw photometry data from multiple recordings into a single
%%structure and parse by specified behavioral state, for FFT analysis with
%%function getFft
%
% [rawS] = extractRaw_fft() - extract new data into structure
% [rawS] = extractRaw_fft(fPath, fName) - extract new data into structure
% [rawS] = extractRaw_fft(rawS) - only parse already extracted data
%
% Description: Extract raw photometry signal from multiple recording files 
% into a larger structure to be used for FFT analysis with function getFft
% User will be prompted to select photometry signal and behavioral state 
% to be used. Recordings lacking behavioral state (e.g. no immobility) 
% will not be included in the output structure.
%
% INPUTS
%   'fPath' - Character array containing folder path where data files are
%       example: 'R:\tritsn01labspace\Anya\FiberPhotometry\AK201-206\220105'
%   'fName' - Cell array, with each cell containing file names for each
%               recording to be added to structure
%   'rawS' - Structure previously generated using extractRaw_fft function
%
% OUPUTS
%   'rawS' - Structure with raw photometry signals from multiple recordings
%       rawS(x).fp_sub, field containing parsed photometry signal during
%       specified behavioral state
%
% Anya Krok, January 2022
%
    %% INPUTS
    switch nargin
        case 0
            [fName,fPath] = uigetfile('*.mat','MultiSelect','On');
            if ~iscell(fName); fName = {fName}; end
        case 2
            fPath = varargin{1};
            fName = varargin{2};
            if ~iscell(fName); fName = {fName}; end
        case 1
            rawS = varargin{1};
    end
    behState = menu('Select behavioral state','Immobility','Locomotion','Reward','Full Trace');
    pickFP = menu('Select photometry signal','ACh','DA');
    
    %% Extract data
    if nargin ~= 1 % If need to load data into new structure
        rawS = struct;
        h = waitbar(0, 'Extracting raw photometry signals into structure');
        for f = 1:length(fName)
            fprintf('Extracting raw photometry data %s ... ',fName{f});
            load(fullfile(fPath,fName{f})); % Load raw data file
            [an,b] = strtok(fName{f},'_'); [day] = strtok(b,'_'); % Parse file name
            x = 1 + length(rawS);
            rawS(x).rec = [an,'-',day]; 
            rawS(x).site = 'DLS'; % CHANGE or remove

            %% Pull parameters required for this analysis
            if isfield(data.gen,'params')
                params = data.gen.params; % Extract params structure
                dsRate = params.dsRate; 
                dsType = params.dsType; % General downsampling parameter
                rawFs = data.gen.acqFs; 
                Fs = data.gen.Fs;
            else
                error('No parameters saved during processData');
            end

            %% Extract photometry and behavior data
            rawS(x).FPnames = data.acq.FPnames;
            rawS(x).rawFP = data.acq.FP;
            rawS(x).rawFs = rawFs;
            if isfield(data,'final')
                if isfield(data.final,'mov')
                    if ~isempty(data.final.mov.onsets)
                        rawS(x).on = data.final.mov.onsets.*dsRate;
                        rawS(x).off = data.final.mov.offsets.*dsRate;
                    end
                    if ~isempty(data.final.mov.onsetsRest)
                        rawS(x).onRest = data.final.mov.onsetsRest.*dsRate;
                        rawS(x).offRest = data.final.mov.offsetsRest.*dsRate;
                    end
                end
            end
            if isfield(data.final,'rew')
                if isfield(data.final.rew,'onset')
                    rawS(x).reward = data.final.rew.onset.*dsRate;
                end
            end
            fprintf('DONE.\n');
            waitbar(f/length(fName),h);

        end
        close(h);
        if isempty(rawS(1).rawFs); rawS(1) = []; end
    end

    %% Parse photometry signal during specified behavioral state
    rmv = zeros(length(rawS),1);
    switch behState
        case 1
            fprintf('Extracting signal during %s (this will take a while!) ...',behState)
            h = waitbar(0,'Extracting signal during behavioral state');
            for x = 1:length(rawS)
                nSampRaw = length(rawS(x).rawFP{pickFP});
                if isfield(rawS,'reward')
                    rewWindow = rawS(x).rawFs;
                    idx_rew = extractEventST([1:nSampRaw]', floor(rawS(x).reward), floor(rawS(x).reward)+rewWindow, 1); % identify recording indices during reward
                else; idx_rew = [];
                end
                idx_imm = extractEventST([1:nSampRaw]', rawS(x).onRest, rawS(x).offRest, 1); % identify recording indices during immobility
                if isempty(idx_imm)
                    fprintf('%s - no immobility \n', rawS(x).rec); 
                    rmv(x) = 1; % change index to be 1, this recording will be removed
                end
                idx_imm = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward from immobility indices
                rawS(x).fp_sub = rawS(x).rawFP{pickFP}(idx_imm); % extract signal during immobility
                rawS(x).fp_lbl = rawS(x).FPnames{pickFP};
                rawS(x).behState = 'immobility';
                waitbar(x/length(rawS),h);
            end; close(h);
        case 2
            fprintf('Extracting signal during %s (this will take a while!) ...',behState)
            h = waitbar(0,'Extracting signal during behavioral state');
            for x = 1:length(rawS)
                nSampRaw = length(rawS(x).rawFP{pickFP});
                if isfield(rawS,'reward')
                    rewWindow = rawS(x).rawFs;
                    idx_rew = extractEventST([1:nSampRaw]', floor(rawS(x).reward), floor(rawS(x).reward)+rewWindow, 1); % identify recording indices during reward
                else; idx_rew = [];
                end
                idx_loc = extractEventST([1:nSampRaw]', rawS(x).on, rawS(x).off, 1); % identify recording indices during locomotion
                if isempty(idx_loc)
                    fprintf('%s - no locomotion \n', rawS(x).rec); 
                    rmv(x) = 1; % change index to be 1, this recording will be removed 
                end
                idx_loc = idx_loc(~ismember(idx_loc, idx_rew)); % exclude reward from locomotion indices
                rawS(x).fp_sub = rawS(x).rawFP{pickFP}(idx_loc); % extract signal during locomotion
                rawS(x).fp_lbl = rawS(x).FPnames{pickFP};
                rawS(x).behState = 'locomotion';
                waitbar(x/length(rawS),h);
            end; close(h);
        case 3
            fprintf('Extracting signal during %s (this will take a while!) ...',behState)
            h = waitbar(0,'Extracting signal during behavioral state');
            for x = 1:length(rawS)
                nSampRaw = length(rawS(x).rawFP{pickFP});
                if isfield(rawS,'reward')
                    rewWindow = rawS(x).rawFs;
                    idx_rew = extractEventST([1:nSampRaw]', floor(rawS(x).reward), floor(rawS(x).reward)+rewWindow, 1); % identify recording indices during reward
                    rawS(x).fp_sub = rawS(x).rawFP{pickFP}(idx_rew); % extract signal during reward
                    rawS(x).fp_lbl = rawS(x).FPnames{pickFP};
                    rawS(x).behState = 'reward';
                else
                    error('no reward');
                end
                waitbar(x/length(rawS),h);
            end; close(h);
        case 4
            for x = 1:length(rawS)
                rawS(x).fp_sub = rawS(x).rawFP{pickFP}; % full trace
                rawS(x).fp_lbl = rawS(x).FPnames{pickFP};
                rawS(x).behState = 'full';
            end
    end
    rawS(rmv == 1) = []; % remove recordings where no fp_sub extracted
    fprintf('DONE! \n');
end