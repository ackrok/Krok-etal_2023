function [mat, staMaxImm, staZeroImm, staShuff50, staShuff95] = getFP2syncSpike(units, beh)
%
% [mat, staMaxImm, staZeroImm, staShuff50, staShuff95] = getFP2syncSpike(units, beh)
%
%
    %% 
    bin = 0.02; window = [-0.01 0.01]; % window for PETH
    uni = unique({units.rec}); % How many unique recordings are there?
    mat = struct; % Initialize structure
    h = waitbar(0, 'Aligning spike times to spike times');
    for u = 1:length(uni)
        ii = find(strcmp({units.rec},uni{u})); % Index of units from this recording
        ib = find(strcmp({beh.rec},uni{u})); % Index of matching photometry data
        if length(ii) < 2; continue; end % Skip recordings where less than 2 units
        sub = units(ii); % New structure with only units from this recording
        for x = 1:length(ii)
            %% Extract spike times
            jj = [1:length(ii)]; jj(x) = []; % "other" 
            st = [sub(x).st]; % Extract spike times of this unit
            st_other = {sub(jj).st}; % Extract spike times of all other units from this recording
            %% PETH
            peth = getClusterPETH(st_other, st, bin, window); % PETH: spike times aligned to spike times
            %%
            prob = []; 
            t_0 = 1; cts_0 = []; cts = {};
            for y = 1:length(st_other)
                cts{y} = peth.cts{y}; 
                prob(:,y) = (sum(peth.cts{y},2))./length(st_other{y}); % Probability of firing = cts/st
                cts_0(y,:) = peth.cts{y}(t_0,:);
            end 
            %% STA
            signal = beh(ib).FP{1}; 
            Fs = beh(ib).Fs;
            signal = signal - nanmean(signal);
            % fp = [beh(1).vel(1); diff(movmean(beh(1).vel,10))]; Fs = 50;
            [sta_fp, sta_time] = getSTA(signal, st, Fs, [-1, 1]);

            %% Load into output structure
            next = 1 + length(mat);
            mat(next).rec = sub(x).rec; 
            mat(next).n = sub(x).n; mat(next).m = [sub(jj).n];
            mat(next).st = st;
            mat(next).cts0 = cts_0; 
            mat(next).cts = cts; mat(next).prob = prob; 
            mat(next).sta = sta_fp;

            %% REST vs MVMT
            stImm = extractEventST(st, beh(ib).onRest, beh(ib).offRest, 1);
            idxImm = find(ismember(st, stImm)); % Idx of spike times during immobility
            stLoc = extractEventST(st, beh(ib).on, beh(ib).off, 1);
            idxLoc = find(ismember(st, stLoc)); % Idx of spike times during locomotion
            mat(next).idxImm = idxImm; 
            mat(next).idxLoc = idxLoc;
            
        end
        waitbar(u/length(uni),h);
    end; close(h);
    if isempty(mat(1).n); mat(1) = []; end
    fprintf('Aligning photometry to spike times complete. \n');
    time = peth.time; 

    %% SHUFF
    nShuff = 10; % change number of shuffles
    for x = 1:length(mat)
        st = mat(x).st;
        stImm = st(mat(x).idxImm); 
        if isempty(stImm); continue; end
        shuffSt = shuffleST(stImm, nShuff); % shuff spike times
        ib = strcmp({beh.rec},mat(x).rec);
        signal = beh(ib).FP{1}; 
        Fs = beh(ib).Fs;
        signal = signal - nanmean(signal);
        prc5 = []; prc50 = []; prc95 = []; % initialize output
        for y = 1:nShuff
            tmp2 = getSTA(signal, shuffSt{y}, Fs, [-1, 1]); % align photometry to shuffled spike times
            prc = prctile(tmp2, [5 50 95], 2);
            prc5(:,y) = prc(:,1);
            prc50(:,y) = prc(:,2);
            prc95(:,y) = prc(:,3);
        end
        tmp = [nanmean(prc5,2), nanmean(prc50,2), nanmean(prc95,2)];
        mat(x).shuffImm = tmp./nShuff; 
    end

    %% (3) 1/N vs N/N IMMOBILITY
    staMax = cell(length(mat),2);
    staMaxImm = nan(size(mat(1).sta,1),length(mat)); 
    staZeroImm = nan(size(mat(1).sta,1),length(mat)); 
    staShuff50 = nan(size(mat(1).sta,1),length(mat)); 
    staShuff95 = nan(size(mat(1).sta,1),length(mat)); 
    propN = [];
    for x = 1:length(mat)
        if isempty(mat(x).idxImm); continue; end
        cts = mat(x).cts0; % extract counts
        cts(cts > 1) = 1; 
        cts = sum(cts, 1); % sum across units
        iN = find(cts >= length(mat(x).m)); iN = iN'; % index of max coherence among units
        i0 = find(cts == 0); i0 = i0'; % index of no coherence
        iN_imm = mat(x).idxImm((ismember(mat(x).idxImm, iN))); % index of max coherence among units that occur during immobility
        i0_imm = mat(x).idxImm((ismember(mat(x).idxImm, i0))); % index of zero coherence among units that occur during immobility
        pull_sta = mat(x).sta;
        pull_sta = pull_sta - nanmean(pull_sta([1:10],:));
        staMax{x,1} = pull_sta(:,i0_imm);
        staMax{x,2} = pull_sta(:,iN_imm);
        staZeroImm(:,x) = nanmean(staMax{x,1},2);
        staMaxImm(:,x) = nanmean(staMax{x,2},2);
        staShuff50(:,x) = mat(x).shuffImm(:,2);
        staShuff95(:,x) = mat(x).shuffImm(:,3) - mat(x).shuffImm(:,2);
        propN(x) = length(iN_imm)/length(cts);
    end
end
