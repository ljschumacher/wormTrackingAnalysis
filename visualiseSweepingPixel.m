% script visualises sweeping using pixel data from selected frames.
% step 1: read frame (or frames, if averaging over several frames); blur frame as specified
% step 2: generate binary image based on intensity threshold to pick out worm pixels from background
% step 3: apply area thresholding to binary image to pick up clusters
% step 4: draw clusters over time

close all
clear

%% set analysis parameters
dataset = 2;
phase = 'fullMovie';
wormnum = '40';
markerType = 'pharynx';
blobAreaThreshold = 8000;
sampleEveryNSec = 120; 
dilateBinaryImage = true;
plotCentroid = false;

numMovieFrames = 3600/sampleEveryNSec;
if dilateBinaryImage
    dilationRadius = 10;
end

%% set fixed parameters

if dataset ==1
    strains = {'npr1','N2','HA'};%{'npr1','HA','N2'}
    assert(~strcmp(markerType,'bodywall'),'Bodywall marker for dataset 1 not available')
elseif dataset ==2
    strains = {'npr1','N2'};
end
pixelsize = 100/19.5; % 100 microns are 19.5 pixels

% filtering parameters
if dataset == 1
    intensityThresholds = containers.Map({'40','HD','1W'},{50, 40, 100});
elseif dataset ==2
    intensityThresholds = containers.Map({'40','HD','1W'},{60, 40, 100});
end
if strcmp(markerType,'pharynx')
    maxBlobSize = 1e5;
    channelStr = 'g';
else
    error('unknown marker type specified, should be pharynx or bodywall')
end

% export fig parameters
exportOptions = struct('Format','EPS2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

%% other initial set up

addpath('auxiliary/')
addpath('visualisation/')

plotColors = parula(numMovieFrames);

%% loop through strains
for strainCtr = 1%:length(strains)
    % load file lists
    if dataset == 1
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_list.xlsx'],1,'A1:E15','basic');
    elseif dataset == 2
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_' channelStr '_list.xlsx'],1,'A1:E15','basic');
    end
    % initialise
    clusterCentroidSpeedFig = figure; hold on
    recordingLabels = cell(1,length(strains));
    recordingColors = distinguishable_colors(length(filenames));
    % loop through files
    for fileCtr = 1:length(filenames) % can be parfor
        filename = filenames{fileCtr};
        recordingLabels{fileCtr} = filename(end-22:end-18);
        % load data
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        trajData = h5read(filename,'/trajectories_data');
        skeletonfilename = strrep(filename,'Results','MaskedVideos');
        skeletonfilename = strrep(skeletonfilename,'_skeletons','');
        fileInfo = h5info(skeletonfilename);
        dims = fileInfo.Datasets(2).Dataspace.Size;
        % get phase restriction frames
        [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr);
        firstFrame = firstFrame+1; % adjust for python zero-indexing
        lastFrame = lastFrame+1;
        % generate a list of movie frame
        movieFrames = round(linspace(firstFrame,double(lastFrame),numMovieFrames));
        % initialise figure and centroid coords/speed matrices
        clusterVisFig = figure; hold on
        blobCentroidCoords = NaN(2,numel(movieFrames));
        clusterCentroidSpeed{fileCtr} = NaN(1,numel(movieFrames)-1);
        % read each movie frame from masked images
        for frameCtr = 1:numel(movieFrames)
            imageFrame = h5read(skeletonfilename,'/mask',[1,1,movieFrames(frameCtr)],[dims(1),dims(2),1]);
            % generate binary segmentation based on black/white contrast
            binaryImage = imageFrame > intensityThresholds(wormnum); % apply blob heat map intensity threshold values
            binaryImage = imfill(binaryImage, 'holes');
            % dilate image
            if dilateBinaryImage
                SE = strel('sphere',dilationRadius);
                binaryImage = imdilate(binaryImage,SE);
            end
            % filter by blob size
            blobMeasurements = regionprops(binaryImage, 'Area','Centroid');
            blobLogInd = [blobMeasurements.Area] > blobAreaThreshold; % apply blob area threshold values
            blobBoundaries = bwboundaries(binaryImage,8,'noholes');
            % plot individual blob boundaries that meet area threshold requirements
            set(0,'CurrentFigure',clusterVisFig)
            for blobCtr = 1:numel(blobLogInd) 
                if blobLogInd(blobCtr)
                    fill(blobBoundaries{blobCtr}(:,1)/size(binaryImage,1)*12,blobBoundaries{blobCtr}(:,2)/size(binaryImage,2)*12,plotColors(frameCtr,:),'edgecolor','none')
                    alpha 0.5
                end
            end
            % get centroids
            if nnz(blobLogInd)>0
                blobCentroidsCoords = reshape([blobMeasurements.Centroid],[2, numel([blobMeasurements.Centroid])/2]);
                %if ~strcmp(filename,'/data2/shared/data/twoColour/Results/recording63/recording63.4g100-250/recording63.4g_X1_skeletons.hdf5')
                % this particular movie has two clusters
                [~,maxAreaIdx] = max([blobMeasurements.Area]);
                blobCentroidCoords(1,frameCtr) = blobCentroidsCoords(1,maxAreaIdx);
                blobCentroidCoords(2,frameCtr) = blobCentroidsCoords(2,maxAreaIdx);
                if plotCentroid
                    set(0,'CurrentFigure',clusterVisFig)
                    plot(blobCentroidCoords(2,frameCtr)/size(binaryImage,1)*12,blobCentroidCoords(1,frameCtr)/size(binaryImage,2)*12,'k--x')
                end
                %end
            end
        end
        % format figure for export
        set(0,'CurrentFigure',clusterVisFig)
        xlim([0 12])
        ylim([0 12])
        xticks([0:2:12])
        yticks([0:2:12])
        xlabel('x (mm)'), ylabel('y (mm)')
        title([strains{strainCtr} ' ' strrep(filename(end-32:end-18),'/','')])
        colorbar
        caxis([0 60])
        cb = colorbar; cb.Label.String = 'minutes';
        if plotCentroid
            figurename = ['figures/sweeping/' strains{strainCtr}...
                '_' strrep(strrep(filename(end-32:end-18),' ',''),'/','') '_blobsOverTimePixelCentroid_blobArea'  num2str(blobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dilateR' num2str(dilationRadius) '_'  phase '_data' num2str(dataset)];
        else
            figurename = ['figures/sweeping/' strains{strainCtr}...
                '_' strrep(strrep(filename(end-32:end-18),' ',''),'/','') '_blobsOverTimePixel_blobArea'  num2str(blobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dilateR' num2str(dilationRadius) '_'  phase '_data' num2str(dataset)];
        end
        exportfig(clusterVisFig,[figurename '.eps'],exportOptions)
        % calculate centroid speed (in microns per minute)
        for frameCtr = 1:numel(movieFrames)-10
            if ~isnan(blobCentroidCoords(1,frameCtr))
                for stepCtr = 1:10 % in case the next sample frame has no cluster, go up to 10 time steps away
                    if ~isnan(blobCentroidCoords(1,frameCtr+stepCtr))
                        break
                    end
                end
                clusterCentroidSpeed{fileCtr}(frameCtr) = sqrt((blobCentroidCoords(1,frameCtr+stepCtr)-blobCentroidCoords(1,frameCtr))^2 +...
                    (blobCentroidCoords(2,frameCtr+stepCtr)-blobCentroidCoords(2,frameCtr))^2)/stepCtr*pixelsize*60/sampleEveryNSec;
            end
        end
        set(0,'CurrentFigure',clusterCentroidSpeedFig)
        recordingsPlotX = 1:sampleEveryNSec/60:60;
        if numel(recordingsPlotX) == numel(clusterCentroidSpeed{fileCtr})
            plot(recordingsPlotX,smoothdata(clusterCentroidSpeed{fileCtr}),'Color',recordingColors(fileCtr,:))
        elseif numel(recordingsPlotX) < numel(clusterCentroidSpeed{fileCtr})
            plot(recordingsPlotX,smoothdata(clusterCentroidSpeed{fileCtr}(1:recordingsPlotX)),'Color',recordingColors(fileCtr,:))
        elseif numel(recordingsPlotX) > numel(clusterCentroidSpeed{fileCtr})
            plot(recordingsPlotX(1:numel(clusterCentroidSpeed{fileCtr})),smoothdata(clusterCentroidSpeed{fileCtr}),'Color',recordingColors(fileCtr,:))
        end
    end
    % export cluster centroid speed figure
    set(0,'CurrentFigure',clusterCentroidSpeedFig)
    xlabel('Time (min)'), ylabel('Cluster Speed (microns/min)')
    legend(recordingLabels,'Location','eastoutside')
    xlim([0 60])
    figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_data' num2str(dataset)];
    exportfig(clusterCentroidSpeedFig,[figurename '.eps'],exportOptions)
    % save cluster centroid speednumbers
    save(['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_data' num2str(dataset) '.mat'],'clusterCentroidSpeed');
    % calculate average speed 
    for fileCtr = 1:length(filenames) %[1,2,3,4,6,10,11] %
        averageSpeed(fileCtr) = nanmedian(clusterCentroidSpeed{fileCtr});
    end
    averageSpeed = median(averageSpeed)
end