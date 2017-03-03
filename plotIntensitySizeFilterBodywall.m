function [ ] =  plotIntensitySizeFilterBodywall(blobFeats,pixelsize,intensityThreshold,maxBlobSize,figurename)
% makes a diagnostic plot of the joint distribution of blob intensity and
% size for purposes of filtering tracked objects
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

SIfig = figure;
histogram2(blobFeats.intensity_mean,blobFeats.area*pixelsize^2,...
    'DisplayStyle','tile','EdgeColor','none','YBinLimits',[0 5e5])
ylabel('area (\mum^2)')
hold on
plot(intensityThreshold*[1 1],[0 maxBlobSize],'r--')
plot([intensityThreshold 255],[maxBlobSize maxBlobSize],'r--')
yyaxis right
histogram(blobFeats.intensity_mean,'DisplayStyle','stairs','Normalization','Probability')
ylabel('P')
xlabel('pixel intensity')
title(figurename,'FontWeight','normal')
set(SIfig,'PaperUnits','centimeters')
filename = ['figures/diagnostics/sizeVintensity_' strrep(strrep(figurename,'/',''),' ','_')];
exportfig(SIfig,[filename '.eps'],exportOptions)
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);
close(SIfig)

end

