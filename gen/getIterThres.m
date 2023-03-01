function thresVec = getIterThres(vel,indVec,winSize,Fs,nStd,flag)
thresVec = zeros(size(indVec));
if flag == 0
    for n = 1:length(indVec)
        ind = indVec(n)-(winSize*Fs);
        if ind<1; ind = 1; end
        thresVec(n) = nStd*std(vel(ind:indVec(n)));
    end
else
    for n = 1:length(indVec)
        ind = indVec(n)+(winSize*Fs);
        if ind<length(vel); ind = length(vel); end
        thresVec(n) = nStd*std(vel(indVec(n):ind));
    end
end
end