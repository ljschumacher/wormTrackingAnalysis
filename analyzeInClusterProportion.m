% plot percentage of worms in cluster, with the option to restrict movies
% to the initial joining phase

clear
close all

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',3);

%% set parameters
dataset = 2;  % '1' or '2'. To specify which dataset to run the script for.
phase = 'fullMovie'; % 'fullMovie', 'joining', or 'sweeping'.
binSeconds = 30; % create bins measured in seconds - only applied to full movie analysis
if dataset ==1
    strains = {'npr1','N2'}; %{'npr1','HA','N2'}
elseif dataset ==2
    strains = {'npr1','N2'}; %{'npr1','N2'}
end
wormnums = {'40'};%{'40','HD'};
saveResults = false;


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
        numFiles = length(filenames);
        
        for fileCtr = 1:numFiles
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            % filter by blob size and intensity
            trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                intensityThresholds(wormnum),maxBlobSize);
            % filter by in-cluster/lone status
            [~, loneWormLogInd, inClusterLogInd,~] = findWormCategory(filename,inClusterNeighbourNum,minNeighbrDist);
            inClusterLogInd(~trajData.filtered)=false;
            loneWormLogInd(~trajData.filtered)=false;
            % filter by worm category
            totalObjInFrame = trajData.frame_number; totalObjInFrame(~trajData.filtered)=NaN;
            inClusterInFrame = trajData.frame_number; inClusterInFrame(~inClusterLogInd)=NaN;
            loneWormInFrame = trajData.frame_number; loneWormInFrame(~loneWormLogInd)=NaN;
            % apply phase restriction
            [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr);
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
            percentInCluster.(strains{strainCtr}){fileCtr} = smoothdata(inClusterPerSecond./totalObjPerSecond*100);
            percentLoneWorm.(strains{strainCtr}){fileCtr} = smoothdata(loneWormPerSecond./totalObjPerSecond*100);
        end
        
        % pool data across movie, plot and save
        percentInCluster.(strains{strainCtr}) = vertcat(percentInCluster.(strains{strainCtr}){:});
        percentInCluster.(strains{strainCtr})(:,1:2) = 0;
        percentLoneWorm.(strains{strainCtr}) = vertcat(percentLoneWorm.(strains{strainCtr}){:});
        percentLoneWorm.(strains{strainCtr})(:,1:2) = 0;
        xtime = ([1:118]/2);
        shadedErrorBar(xtime,median(percentInCluster.(strains{strainCtr}),1),std(percentInCluster.(strains{strainCtr}),1),'k');
        title([strains{strainCtr} '\_' wormnums{numCtr} '\_inCluster'],'FontWeight','normal')
        xlabel('time(min)')
        ylabel('percentage in cluster')
        xlim([0 60])
        ylim([0 100])
        figurename = (['figures/inClusterProportion/inClusterProportion_' strains{strainCtr} '_' wormnums{numCtr} '_' phase '_data' num2str(dataset) '_pool']);
        inClusterProportionFig  = gcf;
        if saveResults
            exportfig(inClusterProportionFig,[figurename '.eps'],exportOptions)
            system(['epstopdf ' figurename '.eps']);
            system(['rm ' figurename '.eps']);
        end
    end
end

% create plot with both strains on
xtime = ([1:118]/2);
figure; hold on
pH(1) = shadedErrorBar(xtime,percentInCluster.(strains{1}),{@mean,@std},'-b',1);
pH(2) = shadedErrorBar(xtime,percentInCluster.(strains{2}),{@mean,@std},'-r',1);
legend([pH(1).mainLine, pH(2).mainLine],strains{1},strains{2});
xlabel('time(min)')
ylabel('percentage in cluster')
xlim([0 60])
ylim([0 100])
figurename = (['figures/inClusterProportion/inClusterProportion_bothStrains_' wormnums{numCtr} '_'  phase '_data' num2str(dataset) '_pool']);
inClusterProportionFig  = gcf;
if saveResults
    exportfig(inClusterProportionFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end

% save data
if saveResults
    filename = ['./figures/inClusterProportion/inClusterProportion_' wormnums{numCtr} '_' phase '_data' num2str(dataset) '.mat'];
    save(filename,'percentInCluster','percentLoneWorm')
end
