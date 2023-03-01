function out = getSynchrony(units, beh)
% Compute pCIN-pCIN cross-correlogram 
%
% out = getSynchrony(units, beh)
%
% INPUTS
%   'units' - structure containing spike times for units and recording ID, 
%       for matching to behavioral data in 'beh'
%   'beh' - structure containing behavioral data, including on and off 
%       times for periods of immobility and locomotion.
%       Must also include recording ID in matching format for matching 
%       behavioral data to spiking data
%
% OUTPUT
%   'out' - structure containing cross-correlograms for all units
%
% Anya Krok, January 2022

%% CCG rest and mvmt spikes
mat = struct; %Initialize structure to save CCG output data into
uni = unique({units.rec}); %Find unique recording IDs across all units
Fs = beh(1).Fs; %Sampling frequency for behavioral data

for x = 1:length(uni)
    fprintf('%s \n',uni{x});
    ii = find(strcmp({units.rec},uni{x})); %Find all units that match this unique recording ID
    iiB = find(strcmp({beh.rec},uni{x})); %Find behavior data that matches this unique recording ID
    if isempty(iiB); continue; end
    if length(ii) < 2; continue; end %If there are <2 units (thus no possible unit pairs), continue to next unique recording ID
    if all(logical(~rem(beh(iiB).onRest,1))); diffFs = 50; else; diffFs = 1; end % adjustment for whether immobility onset/offset times are in seconds or samples
    ccgst = rununitccg({units(ii).st},...
        beh(iiB).on/diffFs, beh(iiB).off/diffFs,...
        beh(iiB).onRest/diffFs, beh(iiB).offRest/diffFs); 

    mat(x).rec = uni{x};
    
    mat(x).n = [ccgst.pairs.n]; mat(x).m = [ccgst.pairs.m];
    mat(x).fr = [];
    mat(x).ccg = [ccgst.pairs.ccg];
    mat(x).shuffPrc = cell(length(ccgst.pairs),1); % percentile of shuffled CCG's - 2.5%, 50%, 97.5%
    
    mat(x).fr_mvmt = [];
    mat(x).ccg_mvmt = [ccgst.pairs.ccg_mvmt]; 
    mat(x).shuffPrc_mvmt = cell(length(ccgst.pairs),1); % percentile of shuffled CCG's - 2.5%, 50%, 97.5%
    
    mat(x).fr_rest = [];
    mat(x).ccg_rest = [ccgst.pairs.ccg_rest]; 
    mat(x).shuffPrc_rest = cell(length(ccgst.pairs),1); % percentile of shuffled CCG's - 2.5%, 50%, 97.5%
    
    mat(x).dist = [];
    for y = 1:length(ccgst.pairs)
        mat(x).fr(y) = 1/mean(diff(ccgst.times(ccgst.pairs(y).m).full));
        mat(x).fr_mvmt(y) = 1/mean(diff(extractEventST(ccgst.times(ccgst.pairs(y).m).full, beh(iiB).on/diffFs, beh(iiB).off/diffFs, 0)));
        mat(x).fr_rest(y) = 1/mean(diff(extractEventST(ccgst.times(ccgst.pairs(y).m).full, beh(iiB).onRest/diffFs, beh(iiB).offRest/diffFs, 0)));
        mat(x).shuffPrc{y} = prctile(ccgst.pairs(y).shuff,[2.5 50 97.5],2); %5th, 50th, 95th percentile of shuffled CCG's
        mat(x).shuffPrc_mvmt{y} = prctile(ccgst.pairs(y).shuff_mvmt,[2.5 50 97.5],2); %5th, 50th, 95th percentile of shuffled CCG's
        mat(x).shuffPrc_rest{y} = prctile(ccgst.pairs(y).shuff_rest,[2.5 50 97.5],2); %5th, 50th, 95th percentile of shuffled CCG's
%         xc_diff = sub(idx(ccgst.pairs(y).n)).coor(1) - sub(idx(ccgst.pairs(y).m)).coor(1); %Distance in x- or y-dimension
%         zc_diff = sub(idx(ccgst.pairs(y).n)).coor(2) - sub(idx(ccgst.pairs(y).m)).coor(2); %Distnace in z-dimension
%         mat(x).dist(y) = hypot(xc_diff,zc_diff);
    end
end
time = ccgst.lag;

%% EXTRACT FROM MAT
time = [-2:0.01:2];
ccg_full = []; ccg_loc = []; ccg_imm = [];
ccgDelta = []; ccg95 = []; ccg50 = [];
ccgDelta_loc = []; ccgDelta_95loc = []; ccgDelta_50loc = [];
ccgDelta_imm = []; ccgDelta_95imm = []; ccgDelta_50imm = [];

for x = 1:length(mat)
    if isempty(mat(x).ccg_rest); continue; end
    ccg_full = [ccg_full, mat(x).ccg];
    ccg_loc = [ccg_loc, mat(x).ccg_mvmt]; 
    ccg_imm = [ccg_imm, mat(x).ccg_rest];
    for y = 1:length(mat(x).fr_rest)
        tmp = (mat(x).ccg(:,y) - mat(x).fr(y))./mat(x).fr(y); %deltaFR = (rate - mu(rate))/mu(rate)
        ccgDelta = [ccgDelta, tmp];
        tmp = (mat(x).ccg_rest(:,y) - mat(x).fr_rest(y))./mat(x).fr_rest(y); %deltaFR = (rate - mu(rate))/mu(rate)
        ccgDelta_imm = [ccgDelta_imm, tmp];
        tmp = (mat(x).ccg_mvmt(:,y) - mat(x).fr_mvmt(y))./mat(x).fr_mvmt(y); %deltaFR = (rate - mu(rate))/mu(rate)
        ccgDelta_loc = [ccgDelta_loc, tmp];
        
        tmp_95 = (mat(x).shuffPrc{y}(:,3) - mat(x).fr(y))./mat(x).fr(y); 
        tmp_50 = (mat(x).shuffPrc{y}(:,2) - mat(x).fr(y))./mat(x).fr(y);
        ccg95 = [ccg95, tmp_95];
        ccg50 = [ccg50, tmp_50];
        
        tmp_95 = (mat(x).shuffPrc_rest{y}(:,3) - mat(x).fr_rest(y))./mat(x).fr_rest(y); 
        tmp_50 = (mat(x).shuffPrc_rest{y}(:,2) - mat(x).fr_rest(y))./mat(x).fr_rest(y);
        ccgDelta_95imm = [ccgDelta_95imm, tmp_95];
        ccgDelta_50imm = [ccgDelta_50imm, tmp_50];
        
        tmp_95 = (mat(x).shuffPrc_mvmt{y}(:,3) - mat(x).fr_mvmt(y))./mat(x).fr_mvmt(y); 
        tmp_50 = (mat(x).shuffPrc_mvmt{y}(:,2) - mat(x).fr_mvmt(y))./mat(x).fr_mvmt(y);
        ccgDelta_95loc = [ccgDelta_95loc, tmp_95];
        ccgDelta_50loc = [ccgDelta_50loc, tmp_50];
    end
end
above95_imm = []; below5_imm = [];
above95_loc = []; below5_loc = [];
for x = 1:length(mat)
    for y = 1:size(mat(x).ccg_rest,2)
        a = mat(x).ccg_rest(:,y); b = mat(x).shuffPrc_rest{y};
        above95_imm = [above95_imm, a > b(:,3)]; %binary vector where CCG passed 95% confidence interval
        below5_imm = [below5_imm, a < b(:,1)]; %binary vector where CCG below 5% confidence interval
    end
    for y = 1:size(mat(x).ccg_mvmt,2)
        a = mat(x).ccg_mvmt(:,y); b = mat(x).shuffPrc_mvmt{y};
        above95_loc = [above95_loc, a > b(:,3)]; %binary vector where CCG passed 95% confidence interval
        below5_loc = [below5_loc, a < b(:,1)]; %binary vector where CCG below 5% confidence interval
    end
end

%% LOAD INTO OUTPUT
out = struct;
y = 1;
out(y).lbl = 'immobility';
out(y).delta = ccgDelta_imm;
out(y).delta50 = ccgDelta_50imm;
out(y).delta95 = ccgDelta_95imm;
out(y).align = ccg_imm;
out(y).above95 = above95_imm;
out(y).below5 = below5_imm;
out(y).time = time;
y = 2;
out(y).lbl = 'locomotion';
out(y).delta = ccgDelta_loc;
out(y).delta50 = ccgDelta_50loc;
out(y).delta95 = ccgDelta_95loc;
out(y).align = ccg_loc;
out(y).above95 = above95_loc;
out(y).below5 = below5_loc;
out(y).time = time;

