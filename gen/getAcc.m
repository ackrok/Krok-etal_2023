function acc = getAcc(vel)
% Simple function for converting velocity signal (cm/s) into acceleration
%
% acc = getAcc(vel)
%
% INPUTS
%   'vel' - vector with velocity signal, usually from data.final.vel or
% 
% OUTPUTS
%   'acc' - vector with acceleration signal
%   
    vel_sm = fliplr(movmean(fliplr(movmean(vel,10)),10)); % Smooth velocity, flip left-right, smooth again, flip back
    acc = [vel(1); diff(vel_sm)]; % Acceleration vector is diff of smoothed velocity vector
end
       