close all
clear

sampleTimeStep = 120; % in seconds
dilateBinaryImage = false;
if dilateBinaryImage
    dilationRadius = 5;
end
blobAreaThreshold = 3000; %5000
plotCentroid = true;

exportOptions = struct('Format','EPS2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

addpath('auxiliary/')
addpath('visualisation/')

[annotationNum,annotationFilenames,~] = xlsread('datalists/BFLongSweeping.xlsx',1,'A1:E40','basic');

for repCtr = 1:6 % go through each recording replicate (6 in total)
    clusterVisFig = figure; hold on
    totalFrames = 0;
    totalSegs = 0;
    for segCtr = 1:5 % go through each hour of the recording replicate (5 hours maximum)
        fileIdx = find(annotationNum(:,1) == repCtr & annotationNum(:,2) == segCtr);
        firstFrame{segCtr} = annotationNum(fileIdx,4)+1; % +1 to adjust for python 0 indexing
        lastFrame{segCtr} = annotationNum(fileIdx,5)+1;
        filename{segCtr} = annotationFilenames{fileIdx};
        if lastFrame{segCtr} - firstFrame{segCtr} > 0 % if this recording has any valid frames
            totalFrames = totalFrames+lastFrame{segCtr}-firstFrame{segCtr}+1;
            totalSegs = totalSegs+1;
        end
    end
    % initialise
    plotColors = parula(ceil(totalFrames/sampleTimeStep/25));
    cumFrame = 0;
    leftoverFrames = 0;
    for segCtr = 1:totalSegs
        skelFilename = strrep(strrep(filename{segCtr},'MaskedVideos','Results'),'.hdf5','_skeletons.hdf5');
        % load data
        frameRate = double(h5readatt(skelFilename,'/plate_worms','expected_fps'));
        trajData = h5read(skelFilename,'/trajectories_data');
        fileInfo = h5info(filename{segCtr});
        dims = fileInfo.Datasets(2).Dataspace.Size;
        if leftoverFrames>0
            assert(firstFrame{segCtr} ==1); % if there are leftover frames from the previous segment, then this segment must start from the very first frame
            firstFrame{segCtr} = sampleTimeStep*frameRate-leftoverFrames+firstFrame{segCtr};
        end
        movieFrames = firstFrame{segCtr}:sampleTimeStep*frameRate:floor((lastFrame{segCtr}-firstFrame{segCtr}+1)/sampleTimeStep/frameRate)*sampleTimeStep*frameRate+firstFrame{segCtr};
        for frameCtr = 1:numel(movieFrames)
            imageFrame = h5read(filename{segCtr},'/mask',[1,1,movieFrames(frameCtr)],[dims(1), dims(2), 1]);
            % generate binary segmentation based on black/white contrast
            binaryImage = imageFrame>0 & imageFrame<100;
            binaryImage = imfill(binaryImage, 'holes');
            % dilate image
            if dilateBinaryImage
                SE = strel('sphere',dilationRadius);
                binaryImage = imdilate(binaryImage,SE);
            end
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
                    fill(blobBoundaries{blobCtr}(:,1),blobBoundaries{blobCtr}(:,2),plotColors(cumFrame+frameCtr,:),'edgecolor','none')
                    alpha 0.5
                end
            end
            % get centroids
            if nnz(blobLogInd)>0
                [~,maxAreaIdx] = max([blobMeasurements.Area].*blobLogInd);
                blobCentroidCoords(1,frameCtr) = blobCentroidsCoords(1,maxAreaIdx);
                blobCentroidCoords(2,frameCtr) = blobCentroidsCoords(2,maxAreaIdx);
                if plotCentroid
                    set(0,'CurrentFigure',clusterVisFig)
                    plot(blobCentroidCoords(2,frameCtr),blobCentroidCoords(1,frameCtr),'k--x')
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
    xlim([0 1600])
    ylim([0 1600])
    xticks([0:400:1600])
    yticks([0:400:1600])
    xlabel('pixels')
    ylabel('pixels')
    if plotCentroid
        figurename = ['figures/sweeping/npr1_clusterCentroidSpeedBF_rep' num2str(repCtr) 'Centroid_sampleTimeStep' num2str(sampleTimeStep) '_blobAreaThreshold' num2str(blobAreaThreshold)];
    else
        figurename = ['figures/sweeping/npr1_clusterCentroidSpeedBF_rep' num2str(repCtr) '_sampleTimeStep' num2str(sampleTimeStep) '_blobAreaThreshold' num2str(blobAreaThreshold)];
    end
    exportfig(clusterVisFig,[figurename '.eps'],exportOptions)
end

