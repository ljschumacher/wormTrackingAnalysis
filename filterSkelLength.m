function [ filtered ] = filterSkelLength(skelData,pixelsize,minSkelLength,...
    maxSkelLength,plotDiagnostics,plotName)
if nargin<6
    plotName = [];
    if nargin<5
        plotDiagnostics = false;
        if nargin<4
            maxSkelLength = Inf;
        end
    end
end

skelLengths = squeeze(sum(sqrt(sum((diff(skelData,1,2)*pixelsize).^2))));
filtered = skelLengths<=maxSkelLength&skelLengths>=minSkelLength;

if plotDiagnostics
    plotSkelLengthDist(skelData(:,:,trajData.filtered),pixelsize,...
        minSkelLength,maxSkelLength,plotName);
end

end

