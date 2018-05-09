close all
clear

sampleEveryNSec = 120; % in seconds
blobAreaThreshold = 3000; % single worm area ~ 500 
plotCentroid = true;

exportOptions = struct('Format','EPS2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

addpath('auxiliary/')

[annotationNum,annotationFilenames,~] = xlsread('datalists/BFLongSweeping.xlsx',1,'A1:E40','basic');

clusterCentroidSpeedFig = figure; hold on
recordingColors = distinguishable_colors(6);

for fileCtr = 1:6 % go through each recording replicate (6 in total)
    clusterVisFig = figure; hold on
    totalFrames = 0;
    totalSegs = 0;
    for segCtr = 1:5 % go through each hour of the recording replicate (5 hours maximum)
        fileIdx = find(annotationNum(:,1) == fileCtr & annotationNum(:,2) == segCtr);
        firstFrame{segCtr} = annotationNum(fileIdx,4)+1; % +1 to adjust for python 0 indexing
        lastFrame{segCtr} = annotationNum(fileIdx,5)+1;
        filename{segCtr} = annotationFilenames{fileIdx};
        if lastFrame{segCtr} - firstFrame{segCtr} > 0 % if this recording has any valid frames
            totalFrames = totalFrames+lastFrame{segCtr}-firstFrame{segCtr}+1;
            totalSegs = totalSegs+1;
        end
    end
    totalSampleFrames = ceil(totalFrames/sampleEveryNSec/25);
    % initialise
    plotColors = parula(totalSampleFrames);
    blobCentroidCoords = NaN(2,totalSampleFrames);
    cumFrame = 0; % keep track of cumulative frames across replicate segments
    leftoverFrames = 0; % keep track of leftover frames at the end of one segment that combines with the start of the next segment
    for segCtr = 1:totalSegs
        skelFilename = strrep(strrep(filename{segCtr},'MaskedVideos','Results'),'.hdf5','_skeletons.hdf5');
        % load data
        frameRate = double(h5readatt(skelFilename,'/plate_worms','expected_fps'));
        pixelsize = double(h5readatt(skelFilename,'/trajectories_data','microns_per_pixel')); % 10 microns per pixel
        trajData = h5read(skelFilename,'/trajectories_data');
        fileInfo = h5info(filename{segCtr});
        dims = fileInfo.Datasets(2).Dataspace.Size;
        if leftoverFrames>0
            assert(firstFrame{segCtr} ==1); % if there are leftover frames from the previous segment, then this segment must start from the very first frame
            firstFrame{segCtr} = sampleEveryNSec*frameRate-leftoverFrames+firstFrame{segCtr};
        end
        movieFrames = firstFrame{segCtr}:sampleEveryNSec*frameRate:...
            floor((lastFrame{segCtr}-firstFrame{segCtr}+1)/sampleEveryNSec/frameRate)*sampleEveryNSec*frameRate+firstFrame{segCtr};
        for frameCtr = 1:numel(movieFrames)
            imageFrame = h5read(filename{segCtr},'/mask',[1,1,movieFrames(frameCtr)],[dims(1),dims(2),1]);
            % generate binary segmentation based on black/white contrast
            binaryImage = imageFrame>0 & imageFrame<100;
            binaryImage = imfill(binaryImage, 'holes');
            % filter by blob size
            blobMeasurements = regionprops(binaryImage, 'Area','Centroid');
            blobCentroidsCoords = reshape([blobMeasurements.Centroid],[2, numel([blobMeasurements.Centroid])/2]);
            blobLogInd = [blobMeasurements.Area] > blobAreaThreshold; % apply blob area threshold values
            blobLogInd = blobLogInd & blobCentroidsCoords(2,:) > 250; % get rid of the annoying box at the edge
            blobBoundaries = bwboundaries(binaryImage,8,'noholes');
            % plot individual blob boundaries that meet area threshold requirements
            set(0,'CurrentFigure',clusterVisFig)
            for blobCtr = 1:numel(blobLogInd)
                if blobLogInd(blobCtr)
                    fill(blobBoundaries{blobCtr}(:,1)*pixelsize/1000,blobBoundaries{blobCtr}(:,2)*pixelsize/1000,plotColors(cumFrame+frameCtr,:),'edgecolor','none')
                    alpha 0.5
                end
            end
            % get centroids
            if nnz(blobLogInd)>0
                [~,maxAreaIdx] = max([blobMeasurements.Area].*blobLogInd);
                blobCentroidCoords(1,cumFrame+frameCtr) = blobCentroidsCoords(1,maxAreaIdx);
                blobCentroidCoords(2,cumFrame+frameCtr) = blobCentroidsCoords(2,maxAreaIdx);
                if plotCentroid
                    set(0,'CurrentFigure',clusterVisFig)
                    plot(blobCentroidCoords(2,cumFrame+frameCtr)*pixelsize/1000,blobCentroidCoords(1,cumFrame+frameCtr)*pixelsize/1000,'k--x')
                end
            end
        end
        cumFrame = cumFrame+numel(movieFrames);
        leftoverFrames = lastFrame{segCtr} - max(movieFrames);
    end
    axis equal
    colorbar
    caxis([0 ceil(totalFrames/25/60)])
    cb = colorbar; cb.Label.String = 'minutes';
    xlim([0 20])
    ylim([0 20])
    xticks([0:5:20])
    yticks([0:5:20])
    xlabel('x (mm)')
    ylabel('y (mm)')
    if plotCentroid
        figurename = ['figures/sweeping/npr1_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixelCentroid_blobArea' num2str(blobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBF'];
    else
        figurename = ['figures/sweeping/npr1_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixel_blobArea' num2str(blobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBF'];
    end
    exportfig(clusterVisFig,[figurename '.eps'],exportOptions)
    % calculate centroid speed (in microns per minute)
    clusterCentroidSpeed{fileCtr} = NaN(1,totalSampleFrames);
    for frameCtr = 1:totalSampleFrames-10
        if ~isnan(blobCentroidCoords(1,frameCtr))
            for stepCtr = 1:10 % in case the next sample frame has no cluster, go up to 10 time steps away
                if ~isnan(blobCentroidCoords(1,frameCtr+stepCtr))
                    break
                end
            end
            clusterCentroidSpeed{fileCtr}(frameCtr) = sqrt((blobCentroidCoords(1,frameCtr+stepCtr)-blobCentroidCoords(1,frameCtr))^2 +...
                (blobCentroidCoords(2,frameCtr+stepCtr)-blobCentroidCoords(2,frameCtr))^2)...
                /stepCtr*pixelsize*60/sampleEveryNSec;
        end
    end
    set(0,'CurrentFigure',clusterCentroidSpeedFig)
    plot(1:sampleEveryNSec/60:ceil(totalFrames/25/60),smoothdata(clusterCentroidSpeed{fileCtr}),'Color',recordingColors(fileCtr,:))
end
% export cluster centroid speed figure
set(0,'CurrentFigure',clusterCentroidSpeedFig)
xlabel('Time (min)'), ylabel('Cluster Speed (microns/min)')
legend({'1','2','3','4','5','6'},'Location','eastoutside')
xlim([0 300])
figurename = 'figures/sweeping/npr1_clusterCentroidSpeed_dataBF';
exportfig(clusterCentroidSpeedFig,[figurename '.eps'],exportOptions)
% save cluster centroid speednumbers
save('figures/sweeping/npr1_clusterCentroidSpeed_dataBF.mat','clusterCentroidSpeed');
% calculate average speed
for fileCtr = [1,2,3,6]%1:6 %[1,2,3,4,6] %
    averageSpeed(fileCtr) = nanmedian(clusterCentroidSpeed{fileCtr});
end
averageSpeed = median(averageSpeed)