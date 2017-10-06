function [] = testHierarchicalClusteringData(dataset,phase,wormnum,linkageMethod)
% calculate speed vs neighbr distance, directional correlation, and
% radial distribution functions
% INPUTS
% dataset: 1 or 2. To specify which dataset to run the script for.
% phase: 'joining', 'fullMovie', or 'sweeping'. Script defines stationary phase as: starts at 10% into the movie, and stops at 60% into the movie (HA and N2) or at specified stopping frames (npr-1).
% wormnum: '40', or 'HD'
% plotDiagnostics: true (default) or false
% OUTPUTS
% none returned, but figures are exported
% issues/to-do:
% - seperate into individual functions for each statistic?

%% set other parameters
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',14,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

if dataset ==1
    strains = {'npr1','N2'};%{'npr1','HA','N2'}
elseif dataset ==2
    strains = {'npr1','N2'};
end

nStrains = length(strains);
if dataset == 1
    intensityThresholds = containers.Map({'40','HD','1W'},{50, 40, 100});
elseif dataset ==2
    intensityThresholds = containers.Map({'40','HD','1W'},{60, 40, 100});
end
maxBlobSize = 1e4;
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
mmConversion = 1e-3;
%% go through strains, densities, movies
clustFig = figure;
subplot(1,2,1), hold on
subplot(1,2,2), hold on
for strainCtr = 1:nStrains
    %% load data
    if dataset == 1
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_list.xlsx'],1,'A1:E15','basic');
    elseif dataset == 2
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_g_list.xlsx'],1,'A1:E15','basic');
    end
    numFiles = length(filenames);
    cophCorr = cell(numFiles,1);
    branchHeights= cell(numFiles,1);
    for fileCtr = 1:numFiles % can be parfor
        filename = filenames{fileCtr};
        trajData = h5read(filename,'/trajectories_data');
        blobFeats = h5read(filename,'/blob_features');
        skelData = h5read(filename,'/skeleton');
        assert(size(skelData,1)==2,['Wrong skeleton size for ' filename])
        assert(size(skelData,2)==2,['Wrong skeleton size for ' filename])
        assert(size(skelData,3)==length(trajData.frame_number),['Wrong number of skeleton frames for ' filename])
        assert(length(blobFeats.velocity_x)==length(trajData.frame_number)&&...
            length(blobFeats.signed_speed)==length(trajData.frame_number),['Wrong number of speed frames for ' filename])
        if all(isnan(skelData(:)))
            warning(['all skeleton are NaN for ' filename])
        end
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        %% randomly sample frames to analyze
        [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr);
        numFrames = round((lastFrame-firstFrame)/frameRate);
        framesAnalyzed = randperm((lastFrame-firstFrame),numFrames) + firstFrame; % randomly sample frames without replacement
        %% filter worms
        trajData.has_skeleton = squeeze(~any(any(isnan(skelData)))); % reset skeleton flag for pharynx data
        trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
            intensityThresholds(wormnum),maxBlobSize)...
            &trajData.has_skeleton;
        % apply phase restriction
        phaseFrameLogInd = trajData.frame_number < lastFrame & trajData.frame_number > firstFrame;
        trajData.filtered(~phaseFrameLogInd)=false;
        %% calculate stats
        cophCorr{fileCtr} = cell(numFrames,1);
        branchHeights{fileCtr}= cell(numFrames,1);
        for frameCtr = 1:numFrames % one may be able to vectorise this
            frame = framesAnalyzed(frameCtr);
            [x ,y] = getWormPositions(trajData, frame, true);
            N = length(x);
            if N>1 % need at least two worms in frame
                pairDists = pdist([x y]).*pixelsize.*mmConversion; % distance between all pairs, in micrometer
                clustTree = linkage(pairDists,linkageMethod);
                cophCorr{fileCtr}{frameCtr} = cophenet(clustTree,pairDists);
                branchHeights{fileCtr}{frameCtr} = clustTree(:,3);
            end
        end
        % pool data from frames
        cophCorr{fileCtr} = vertcat(cophCorr{fileCtr}{:});
        branchHeights{fileCtr} = vertcat(branchHeights{fileCtr}{:});
    end
    %% combine data from multiple files
    cophCorr = vertcat(cophCorr{:});
    branchHeights = vertcat(branchHeights{:});
    %% plot histogram of branch heights
    subplot(1,2,1)
    histogram(branchHeights,'Normalization','probability','EdgeColor','none','BinWidth',0.05)
    subplot(1,2,2)
    [hcounts, hedges] = histcounts(cophCorr,'Normalization','probability','BinWidth',0.02);
    stairs(hedges,[1 cumsum(hcounts,'reverse')])
end
%% format and export figures
set(clustFig,'PaperUnits','centimeters')

subplot(1,2,1)
box on
xlim([0 8])
ylim([0 0.1])
legend(strains)
ylabel('P')
xlabel('inter-cluster distance')
subplot(1,2,2)
box on
xlabel('cophenetic correlation')
xlim([0.5 1])
title([linkageMethod ' linkage'],'FontWeight','normal')
figurename = ['figures/clustering/branchHeights_' linkageMethod '_' wormnum '_' phase '_data' num2str(dataset) '_jointraj'];
exportfig(clustFig,[figurename '.eps'],exportOptions)
system(['epstopdf ' figurename '.eps']);
system(['rm ' figurename '.eps']);
%
