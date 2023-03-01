function data = data_dsRateDivisible(data)
%Validate whether signal is divisible by downsampling rate, will be false
%if wavesurfer crashed during recording
%
%   data = data_dsRateDivisible(data)
% 
%   Anya Krok, September 2022
%

    thisLength = length(data.acq.time); % Length of signal that was acquired
    dsRate = data.gen.params.dsRate; % Downsample rate, will only be present if have already run processData
    if floor(thisLength/dsRate) ~= thisLength/dsRate % Check if length is divisible by downsample rate
        setLength = floor(thisLength/dsRate)*dsRate; % Adjust length of all signals to be divisible by dsRate
        data.acq.time = data.acq.time(1:setLength);
        if isfield(data.acq,'wheel')
            data.acq.wheel = data.acq.wheel(1:setLength); end
        if isfield(data.acq,'lick')
            data.acq.lick{1} = data.acq.lick{1}(1:setLength); end
        if isfield(data.acq,'rew')
            data.acq.rew{1} = data.acq.rew{1}(1:setLength); end
        if isfield(data.acq,'FP')
            for y = 1:length(data.acq.FP)
                data.acq.FP{y} = data.acq.FP{y}(1:setLength); end; end
        if isfield(data.acq,'refSig')
            for y = 1:length(data.acq.FP)
                data.acq.refSig{y} = data.acq.refSig{y}(1:setLength); end; end
    end

end