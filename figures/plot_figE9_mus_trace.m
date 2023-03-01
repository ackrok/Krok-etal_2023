%AK236 muscimol

% demodC = cell(1,2);
% for y = 1:2
%     sigEdge = data.gen.params.FP.sigEdge; 
%     rawFs = 5000;
%     rawFP = data.acq.FP{y}; %Extract FP trace
%     modFreq = data.gen.params.FP.modFreq(y); % ANYA EDIT 21/04/06
%     ref = findRef(modFreq,data.acq.refSig,rawFs); %Find the reference signal from the refsig array using modulation frequency
%     demod = digitalLIA(rawFP,ref,modFreq,rawFs,10,8); %Peform the demodulation
%     demod = demod((sigEdge*rawFs)+1:end-(sigEdge*rawFs));
%     demodC{y} = demod;
% end
% fprintf('Done\n');

load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_snc_mus_demod');

%%
for y = 1:2
demod = demodC{y};
clr = {'g','m'};

rawFs = 5000;
fig = figure; fig.Position(3) = 1375;
subplot(1,3,1);
r = [1/rawFs 1200].*rawFs;
plot([r(1)/rawFs:50/rawFs:r(2)/rawFs],demod(r(1):50:r(2)), clr{y});
xlabel('Time (s)'); ylabel('a.u.'); axis square
title(sprintf('%s [0 %d]s',lbl{y},r(2)/rawFs));
switch y; case 2; ylim([0.2 0.4]); 
    case 1; ylim([0.1 0.3]); end

subplot(1,3,2);
r = [1 21]; r = r.*rawFs;
plot([r(1)/rawFs:1/rawFs:r(2)/rawFs],demod(r(1):r(2)), clr{y});
xlabel('Time (s)'); ylabel('a.u.'); axis square
title(sprintf('[%d %d]s',r(1)/rawFs,r(2)/rawFs));
xlim([1 21]);
switch y; case 2; ylim([0.28 0.38]); 
    case 1; ylim([0.15 0.25]); end

subplot(1,3,3);
r = [740 760]; r = r.*rawFs;
% r = [880 900]; r = r.*rawFs;
plot([r(1)/rawFs:1/rawFs:r(2)/rawFs],demod(r(1):r(2)), clr{y});
xlabel('Time (s)'); ylabel('a.u.'); axis square
title(sprintf('[%d %d]s',r(1)/rawFs,r(2)/rawFs));
switch y; case 2; ylim([0.2 0.3]); 
    case 1; ylim([0.16 0.26]); end

movegui(gcf,'center');
end