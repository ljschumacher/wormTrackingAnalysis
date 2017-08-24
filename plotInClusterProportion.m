% plot percentage of worms in cluster, with the option to restrict movies
% to the initial joining phase

clear
close all

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

%% set parameters
dataset = 1;  % '1' or '2'. To specify which dataset to run the script for.
phase = 'sweeping'; % 'fullMovie', 'joining', or 'sweeping'.
binSeconds = 30; % create bins measured in seconds - only applied to full movie analysis
if dataset ==1
    strains = {'N2','npr1'}; %{'npr1','HA','N2'}
elseif dataset ==2
    strains = {'N2','npr1'}; %{'npr1','N2'}
end
wormnums = {'40'};%{'40','HD'};
if dataset == 1
    intensityThresholds = containers.Map({'40','HD','1W'},{50, 40, 100});
elseif dataset ==2
    intensityThresholds = containers.Map({'40','HD','1W'},{60, 40, 100});
end
maxBlobSize = 1e4;
minNeighbrDist = 2000;
inClusterNeighbourNum = 3;
pixelsize = 100/19.5; % 100 microns are 19.5 pixels


%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        % load data
        if dataset == 1
            [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_list.xlsx'],1,'A1:E15','basic');
        elseif dataset == 2
            [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_g_list.xlsx'],1,'A1:E15','basic');
        end
        phaseFrames = phaseFrames-1; % to correct for python indexing at 0
        numFiles = length(filenames);
        inClusterProportionFig = figure; hold on
        loneWormProportionFig = figure; hold on
        for fileCtr = 1:numFiles
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            if strcmp(phase, 'fullMovie')
                firstFrame = 0;
                lastFrame = phaseFrames(fileCtr,4);
            elseif strcmp(phase,'joining')
                firstFrame = phaseFrames(fileCtr,1);
                lastFrame = phaseFrames(fileCtr,2);
            elseif strcmp(phase,'sweeping')
                firstFrame = phaseFrames(fileCtr,3);
                lastFrame = phaseFrames(fileCtr,4);
            end
            % filter by blob size and intensity
            trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                intensityThresholds(wormnum),maxBlobSize);
            % filter by in-cluster/lone status
            min_neighbr_dist = h5read(filename,'/min_neighbr_dist');
            num_close_neighbrs = h5read(filename,'/num_close_neighbrs');
            neighbr_dist = h5read(filename,'/neighbr_distances');
            inClusterLogInd = num_close_neighbrs>=inClusterNeighbourNum;
            inClusterLogInd(~trajData.filtered)=false;
            loneWormLogInd = min_neighbr_dist>=minNeighbrDist;
            loneWormLogInd(~trajData.filtered)=false;
            % filter by worm category
            totalObjInFrame = trajData.frame_number; totalObjInFrame(~trajData.filtered)=NaN;
            inClusterInFrame = trajData.frame_number; inClusterInFrame(~inClusterLogInd)=NaN;
            loneWormInFrame = trajData.frame_number; loneWormInFrame(~loneWormLogInd)=NaN;
            % apply phase restriction
            % firstFrame = double(round(max(trajData.frame_number)/10)); % cut out the first 10 percent of the movie for joining phase restriction
            phaseFrameLogInd = trajData.frame_number <= lastFrame & trajData.frame_number >= firstFrame;
            totalObjInFrame(~phaseFrameLogInd) = false;
            inClusterInFrame(~phaseFrameLogInd) = false;
            loneWormInFrame(~phaseFrameLogInd) = false;
            % binning and extracting counts for each bin
            if strcmp(phase,'fullMovie')
                numBins = (lastFrame-firstFrame)/frameRate/binSeconds;
            else
                numBins = 120; % create 120 bins over restricted phases
            end
            binEdges = linspace(firstFrame,lastFrame,numBins);
            [totalObjPerSecond,~]=histcounts(totalObjInFrame,binEdges);
            [inClusterPerSecond,~]=histcounts(inClusterInFrame,binEdges);
            [loneWormPerSecond,~]=histcounts(loneWormInFrame,binEdges);
            % plot in-cluster/lone worms as percentage of total worms
            % across the movie
            percentInCluster = smooth(inClusterPerSecond./totalObjPerSecond*100,'moving');
            percentLoneWorm = smooth(loneWormPerSecond./totalObjPerSecond*100,'moving');
            set(0,'CurrentFigure',inClusterProportionFig)
            plot(percentInCluster)
            set(0,'CurrentFigure',loneWormProportionFig)
            plot(percentLoneWorm)
        end
        %format plot and export
        set(0,'CurrentFigure',inClusterProportionFig)
        title([strains{strainCtr} '\_' wormnums{numCtr} '\_inCluster'],'FontWeight','normal')
        xlabel('time')
        if strcmp(phase,'fullMovie')
            xlim([0 (3600/binSeconds)]) % set x axis for 1 hour
        end
        ylabel('percentage')
        ylim([0 100])
        set(inClusterProportionFig,'PaperUnits','centimeters')
        figurename = ['figures/inClusterProportion/inClusterProportion_' strains{strainCtr} '_' wormnums{numCtr} '_' phase '_data' num2str(dataset)];
        exportfig(inClusterProportionFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        %
        set(0,'CurrentFigure',loneWormProportionFig)
        title([strains{strainCtr} '\_' wormnums{numCtr} '\_loneWorms'],'FontWeight','normal')
        xlabel('time')
        if strcmp(phase, 'fullMovie')
            xlim([0 (3600/binSeconds)]) % set x axis for 1 hour
        end
        ylabel('percentage')
        ylim([0 50])
        set(loneWormProportionFig,'PaperUnits','centimeters')
        figurename = ['figures/inClusterProportion/loneWormsProportion_' strains{strainCtr} '_' wormnums{numCtr} '_' phase '_data' num2str(dataset)];
        exportfig(loneWormProportionFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
    end
end