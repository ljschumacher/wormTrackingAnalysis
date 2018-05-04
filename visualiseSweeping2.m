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
blobAreaThreshold = 6000; %800;
sampleEveryNSec = 30; 
blurImage = false;
poolHalfSecondFrames = false;
dilateBinaryImage = true;


numMovieFrames = 3600/sampleEveryNSec;
if blurImage
    blurAlpha = 8;
elseif dilateBinaryImage
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
for strainCtr = 1:length(strains)
    % load file lists
    if dataset == 1
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_list.xlsx'],1,'A1:E15','basic');
    elseif dataset == 2
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_' channelStr '_list.xlsx'],1,'A1:E15','basic');
    end
    % loop through files
    for fileCtr = 1:length(filenames) % can be parfor
        filename = filenames{fileCtr};
        % load data
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        if poolHalfSecondFrames
            poolFrameNum = ceil(frameRate/(60/sampleEveryNSec));
        end
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
        clusterVisFig = figure; hold on
        % read each movie frame from masked images
        for frameCtr = 1:numel(movieFrames)
            imageFrame = h5read(skeletonfilename,'/mask',[1,1,movieFrames(frameCtr)],[dims(1), dims(2), 1]);
            if poolHalfSecondFrames
                for poolFrameCtr = 1:poolFrameNum-1
                    if movieFrames(frameCtr)+poolFrameCtr < lastFrame
                        newImageFrame = h5read(skeletonfilename,'/mask',[1,1,movieFrames(frameCtr)+poolFrameCtr],[dims(1), dims(2), 1]);
                    else
                        newImageFrame = h5read(skeletonfilename,'/mask',[1,1,lastFrame],[dims(1), dims(2), 1]);
                    end
                    imageFrame = imageFrame+newImageFrame;
                end
            end
            if blurImage
                imageFrame = imgaussfilt(imageFrame,blurAlpha);
                if frameRate == 9
                    binaryImage = imageFrame >200;
                elseif frameRate == 3
                    binaryImage = imageFrame >100;
                end
            else
                binaryImage = imageFrame > intensityThresholds(wormnum); % apply blob heat map intensity threshold values
            end
            % generate binary segmentation based on black/white contrast
            binaryImage = imfill(binaryImage, 'holes');
            % dilate image
            if dilateBinaryImage
                SE = strel('sphere',dilationRadius);
                binaryImage = imdilate(binaryImage,SE);
            end
            % filter by blob size
            blobMeasurements = regionprops(binaryImage, 'Area');
            blobLogInd = [blobMeasurements.Area] > blobAreaThreshold; % apply blob area threshold values
            blobBoundaries = bwboundaries(binaryImage,8,'noholes');
            [frameCtr nnz(blobLogInd)]
            set(0,'CurrentFigure',clusterVisFig)
            for blobCtr = 1:numel(blobLogInd) % plot individual blob boundaries that meet area threshold requirements
                if blobLogInd(blobCtr)
                    fill(blobBoundaries{blobCtr}(:,1)/size(binaryImage,1)*12,blobBoundaries{blobCtr}(:,2)/size(binaryImage,2)*12,plotColors(frameCtr,:),'edgecolor','none')
                    alpha 0.5
                end
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
        figurename = ['figures/sweeping/' strains{strainCtr}...
            '_' strrep(strrep(filename(end-32:end-18),' ',''),'/','') '_blobsOverTimePixel_blobArea'  num2str(blobAreaThreshold) '_totalFrames' num2str(numMovieFrames) '_dilateR' num2str(dilationRadius) '_'  phase '_data' num2str(dataset)];
        %exportfig(clusterVisFig,[figurename '.eps'],exportOptions)
    end
end