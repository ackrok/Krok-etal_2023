function idxStates = extractBehavioralStates(beh)
%Extract behavioral states using data in beh structure. 
%
% First column contains indices of samples during immobility, excluding 
%   reward if present.
% Second column contains indices of samples during locomotion, also
%   excluding reward if present.
% Third column contains indices of samples during reward delivery period,
%   for rewarded trials, if reward delivery occured.
%
% idxStates = extractBehavioralStates(beh)
%
% Anya Krok, September 2022
%
    idxStates = cell(length(beh),3);
    for x = 1:length(beh)
        idxImm = []; idxMov = []; idxRew = [];
        nSamp = length(beh(x).time); % number of total samples
        Fs = beh(x).Fs; % sampling rate
        if isfield(beh,'on')
            if ~isempty(beh(x).on)
                idxMov = extractEventST([1:nSamp]', beh(x).on, beh(x).off, 1); % identify samples during locomotion
            end
        end
        if isfield(beh,'onRest')
            if ~isempty(beh(x).onRest)
                idxImm = extractEventST([1:nSamp]', beh(x).onRest, beh(x).offRest, 1); % identify samples during locomotion
            end
        end
        if isfield(beh,'reward')
            if ~isempty(beh(x).reward)
                idxRew = extractEventST([1:nSamp]', floor(beh(x).reward), floor(beh(x).reward)+(2*Fs), 1); % identify sample during reward window: [0 +1]seconds relative to reward delivery for all trials (rewarded or not)
                idxMov = idxMov(~ismember(idxMov, idxRew)); % exclude all reward windows from locomotion
                idxImm = idxImm(~ismember(idxImm, idxRew)); % exclude all reward windows from immobility
                rewYes = extractRewardedTrials(beh(x).reward./Fs, beh(x).lick./Fs); % identify rewarded trials
                idxRew = extractEventST([1:nSamp]', floor(beh(x).reward(rewYes))-Fs, floor(beh(x).reward(rewYes))+(2*Fs), 1); % identify sample during reward window for only rewarded trials
            end
        end
        idxStates{x,1} = idxImm; % load into output structure
        idxStates{x,2} = idxMov; % load into output structure 
        idxStates{x,3} = idxRew; % load into output structure
        for z = 1:3
            needL = 200.*Fs;
            idxTmp = idxStates{x,z};
            if ~isempty(idxTmp)
                if length(idxTmp) < needL % ensure that have at least 200 seconds of signal
                    idxTmp = repmat(idxTmp, [ceil(needL/length(idxTmp)) 1]); % duplicate indices to lengthen signal for processing
                    idxStates{x,z} = idxTmp; % re-insert into output structure
                end
            end
        end
    end