function [ filtered ] = filterIntensityAndSize(blobFeats,pixelsize,intensityThreshold,...
    maxBlobSize,plotDiagnostics,plotName)
if nargin<6
    plotName = [];
    if nargin<5
        plotDiagnostics = false;
        if nargin<4
            maxBlobSize = Inf;
        end
    end
end

filtered = (blobFeats.area*pixelsize^2<=maxBlobSize)&...
                    (blobFeats.intensity_mean>=intensityThreshold);

if plotDiagnostics
    plotIntensitySizeFilter(blobFeats,pixelsize,...
        intensityThreshold,maxBlobSize,plotName)
end

end

