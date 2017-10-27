clear
close all

%% set parameters
phase = 'joining'; % 'fullMovie', 'joining', or 'sweeping'.
dataset = 2; % 1 or 2
<<<<<<< HEAD
marker = 'bodywall'; % 'pharynx' or 'bodywall'
=======
marker = 'pharynx'; % 'pharynx' or 'bodywall'
>>>>>>> 4ab6b71013bf87820e87a317cbde0cbe6dab5a8c
strains = {'npr1'}; % {'npr1','N2'}
wormnums = {'40'};% {'40'};
postExitDuration = 5; % duration (in seconds) after a worm exits a cluster to be included in the leave cluster analysis
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
<<<<<<< HEAD
saveResults = true;
useManualTraj = true; % option to use manually joined trajectories; only can be true if using bodywall data
if useManualTraj
    assert(strcmp(marker,'bodywall'),'useManualTraj must only be used in conjunction with bodywall marker')
end
=======
saveResults = false;
>>>>>>> 4ab6b71013bf87820e87a317cbde0cbe6dab5a8c

if dataset == 1
    intensityThresholds_g = containers.Map({'40','HD','1W'},{50, 40, 100});
elseif dataset ==2
    intensityThresholds_g = containers.Map({'40','HD','1W'},{60, 40, 100});
end
intensityThresholds_r = containers.Map({'40','HD','1W'},{60, 40, 100});
maxBlobSize_g = 1e4;
maxBlobSize_r = 2.5e5;
minSkelLength_r = 850;
maxSkelLength_r = 1500;
minNeighbrDist = 2000;
inClusterNeighbourNum = 3;
pixelsize = 100/19.5; % 100 microns are 19.5 pixels

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        % load file list
        if dataset ==1 & strcmp(marker,'pharynx')
            [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_list.xlsx'],1,'A1:E15','basic');
        elseif dataset ==2 & strcmp(marker,'pharynx')
            [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_g_list.xlsx'],1,'A1:E15','basic');
        elseif dataset ==2 & strcmp(marker,'bodywall')
            [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_r_list.xlsx'],1,'A1:E15','basic');
        else
            warning('specified dataset/marker combination does not exist')
        end
        numFiles = length(filenames);
        leaveClusterFilterEffect = NaN(numFiles,6);
        
        %% go through individual movies
<<<<<<< HEAD
        for fileCtr = 1:numFiles
=======
        for fileCtr = 2:numFiles
>>>>>>> 4ab6b71013bf87820e87a317cbde0cbe6dab5a8c
            %% load data
            filename = filenames{fileCtr}
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton'); % in pixels
            
            %% filter worms by various criteria
            if strcmp(marker, 'pharynx')
                % filter green by blob size and intensity
                trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                    intensityThresholds_g(wormnum),maxBlobSize_g);
            elseif strcmp(marker, 'bodywall')
<<<<<<< HEAD
                % filter red by manually joined traj
                if useManualTraj
                    features = h5read(strrep(filename,'skeletons','feat_manual'),'/features_timeseries');
                    trajData.filtered = ismember(trajData.worm_index_manual,int32(features.worm_index));
                end
=======
>>>>>>> 4ab6b71013bf87820e87a317cbde0cbe6dab5a8c
                % filter red by blob size and intensity
                if contains(filename,'55')||contains(filename,'54')
                    intensityThreshold = 80;
                else
                    intensityThreshold = 40;
                end
                trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                    intensityThreshold,maxBlobSize_r);
                % filter red by skeleton length
                trajData.filtered = trajData.filtered&logical(trajData.is_good_skel)...
                    &filterSkelLength(skelData,pixelsize,minSkelLength_r,maxSkelLength_r);
            end
            % apply phase restriction
            [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr);
            phaseFrameLogInd = trajData.frame_number <= lastFrame & trajData.frame_number >= firstFrame;
            trajData.filtered(~phaseFrameLogInd) = false;
            % restrict to only forward-moving worms (bodywall data doesn't currently have signed_speed field)
            if strcmp(marker,'pharynx')
                signedSpeedLogInd = blobFeats.signed_speed>=0;
                trajData.filtered(~signedSpeedLogInd) = false;
            end
            % find worms that have just left a cluster vs lone worms
            trajDataWormIndexJoined = trajData.worm_index_joined;
            trajDataFiltered = trajData.filtered;
            leaveClusterFilterEffect(fileCtr,:) = findWormCategoryFilterEffect(filename,inClusterNeighbourNum,minNeighbrDist,postExitDuration,trajDataWormIndexJoined,trajDataFiltered);
        end
        if saveResults == true
<<<<<<< HEAD
            if useManualTraj
                trajFileName = 'manualTraj_';
            else
                trajFileName = '';
            end
            save(['figures/turns/results/leaveClusterFilterEffect_' trajFileName strain '_' wormnum '_' phase '_data' num2str(dataset) '_' marker '.mat'],'leaveClusterFilterEffect')
=======
            save(['figures/turns/results/leaveClusterFilterEffect_' strain '_' wormnum '_' phase '_data' num2str(dataset) '_' marker '.mat'],'leaveClusterFilterEffect')
>>>>>>> 4ab6b71013bf87820e87a317cbde0cbe6dab5a8c
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local function %%%%%%%%%%%%%%

function leaveClusterFilterEffect = findWormCategoryFilterEffect(filename,inClusterNeighbourNum,minNeighbrDist,postExitDuration,trajDataWormIndexJoined,trajDataFiltered)

% function calculates % of worm frames and worm trajectories thrown out
% during various steps of filtering for leave cluster frames

%% INPUTS:
% filename: full path to the _skeletons.hdf5 file
% inClusterNeighbourNum: [1x1] double = 3. The number of close neighbors needed for worm to be considered 'in cluster'
% minNeighbrDist: [1x1] double = 2000. The minimum distance required for a worm to be considered 'lone worm'
% postExitDuration: [1x1] double = 5. The duration (in seconds) after a
% worm exits a cluster to be included in the leave cluster analysis. Input
% only needed if needing leaveClusterLogInd.

%% OUTPUTS:
% leaveClusterFilterEffect: 6-element row vector. column 1: initial total frame
% count; column 2: percentage of frames filtered by cluster re-entrance;
% column 3: percentage of frames filtered by leaving cluster for good; 4:
% initial total worm count; column 5:  percentage of worms filtered by cluster re-entrance;
% column 6: percentage of worms filtered by leaving cluster for good.

%% load file
min_neighbr_dist = h5read(filename,'/min_neighbr_dist');
num_close_neighbrs = h5read(filename,'/num_close_neighbrs');
neighbr_dist = h5read(filename,'/neighbr_distances');
frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
trajData = h5read(filename,'/trajectories_data');

%% classify worms
%% identify in cluster worms
inClusterLogInd = num_close_neighbrs>=inClusterNeighbourNum;
%% identify lone worms
loneWormLogInd = min_neighbr_dist>=minNeighbrDist;
%% identify leave cluster worms
% find worm-frames where inCluster changes from true to false
leaveClusterLogInd = vertcat(false,inClusterLogInd(1:end-1)&~inClusterLogInd(2:end));
leaveClusterStart = find(leaveClusterLogInd);
% loop through each exit event, retain frames for the specified duration after a worm exits cluster
for exitCtr = 1:numel(leaveClusterStart)
    thisExitIdx = leaveClusterStart(exitCtr);
    wormIndex = trajData.worm_index_joined(thisExitIdx);
    % check for the number of frames that the same worm has beyond the point of cluster exit
    wormPathLength = nnz(trajData.worm_index_joined(thisExitIdx:end)==wormIndex);
    if wormPathLength>=postExitDuration*frameRate
        leaveClusterEnd = leaveClusterStart+postExitDuration*frameRate-1;
    else
        leaveClusterEnd = leaveClusterStart+wormPathLength-1;
    end % this also excludes movie segments with ending frames beyond highest frame number
    % go through each starting frame to generate logical index for leave cluster worms
    leaveClusterLogInd(leaveClusterStart(exitCtr):leaveClusterEnd(exitCtr))=true;
end
initialFrameCount = nnz(leaveClusterLogInd & trajDataFiltered); %
initialWormCount = numel(unique(trajDataWormIndexJoined(trajDataFiltered & leaveClusterLogInd))); %
% exclude when worms move back into a cluster
leaveClusterLogInd(inClusterLogInd)=false;
reEnterFrameCount = initialFrameCount-nnz(leaveClusterLogInd & trajDataFiltered); %
reEnterFramesPercent = reEnterFrameCount/initialFrameCount*100; %
reEnterWormCount = initialWormCount - numel(unique(trajDataWormIndexJoined(trajDataFiltered & leaveClusterLogInd))); %
reEnterWormPercent = reEnterWormCount/initialWormCount * 100; %
% exclude worms that have become lone worm
leaveClusterLogInd(loneWormLogInd)=false;
leaveFrameCount = initialFrameCount-reEnterFrameCount-nnz(leaveClusterLogInd & trajDataFiltered); %
leaveFramesPercent = leaveFrameCount/initialFrameCount*100; %
leaveWormCount = initialWormCount - reEnterWormCount - numel(unique(trajDataWormIndexJoined(trajDataFiltered & leaveClusterLogInd))); %
leaveWormPercent = leaveWormCount/initialWormCount * 100; %
% create output
leaveClusterFilterEffect = [initialFrameCount,reEnterFramesPercent,leaveFramesPercent,initialWormCount,reEnterWormPercent,leaveWormPercent]; %
end