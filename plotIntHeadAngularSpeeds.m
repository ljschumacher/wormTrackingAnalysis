% plot two measures of head angular speed based on red worm data from the
% second dataset: integrated vs. cumulative HAS, to see whether either
% measure is useful for distinguishing between different strains or worm
% category

% issues to consider: check the way consecutive 5s post-exit is selected
% for lone worms - removing NaN breaks up the consecutiveness
% need to add cumulative measure

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
phase = 'sweeping'; % 'fullMovie', 'joining', or 'sweeping'.
strains = {'npr1','N2'}; % {'npr1','N2'}
wormnums = {'40'};% {'40','HD'};
postExitDuration = 5; % set the duration (in seconds) after a worm exits a cluster to be included in the analysis

intensityThresholds = containers.Map({'40','HD','1W'},{60, 40, 100});
maxBlobSize = 2.5e5;
minSkelLength = 850;
maxSkelLength = 1500;
minNeighbrDist = 2000;
inClusterNeighbourNum = 3;
pixelsize = 100/19.5; % 100 microns are 19.5 pixels

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        % load data
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_r_list.xlsx'],1,'A1:E15','basic');
        phaseFrames = phaseFrames-1; % to correct for python indexing at 0
        numFiles = length(filenames);
        % create cell arrays to hold individual movie values to be pooled
        iHAS_leaveCluster_poolFiles = cell(1,numFiles); 
        iHAS_loneWorm_poolFiles = cell(1,numFiles); 
        iHASFig = figure; hold on % iHAS = integrated Head Angular Speed
        %cHASFig = figure; hold on % cHAS = cumulative Head Angular Speed
        %% go through individual movies
        for fileCtr = 1:numFiles
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton');
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            features = h5read(strrep(filename,'skeletons','features'),'/features_timeseries');
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
            %% filter worms by various criteria
            % filter red by blob size and intensity
            if contains(filename,'55')||contains(filename,'54')
                intensityThreshold = 80;
            else
                intensityThreshold = 40;
            end
            trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                intensityThreshold,maxBlobSize);
            % filter red by skeleton length
            trajData.filtered = trajData.filtered&logical(trajData.is_good_skel)...
                &filterSkelLength(skelData,pixelsize,minSkelLength,maxSkelLength);
            % apply phase restriction
            phaseFrameLogInd = trajData.frame_number <= lastFrame & trajData.frame_number >= firstFrame;
            trajData.filtered(~phaseFrameLogInd) = false;
            features.filtered = ismember(features.skeleton_id+1,find(trajData.filtered)); % use trajData.filtered to filter out unwa
            % find worms that have just left a cluster
            min_neighbr_dist = h5read(filename,'/min_neighbr_dist');
            num_close_neighbrs = h5read(filename,'/num_close_neighbrs');
            neighbr_dist = h5read(filename,'/neighbr_distances');
            inClusterLogInd = num_close_neighbrs>=inClusterNeighbourNum;
            leaveClusterLogInd = vertcat(false,inClusterLogInd(1:end-1)&~inClusterLogInd(2:end)); % find worm-frames where inCluster changes from true to false
            leaveClusterFrameStart = find(leaveClusterLogInd);
            leaveClusterFrameEnd = leaveClusterFrameStart+postExitDuration*frameRate;
            leaveClusterFrameEnd = leaveClusterFrameEnd(leaveClusterFrameEnd<=numel(leaveClusterLogInd)); % exclude movie segments with frames beyond highest frame number
            leaveClusterFrameStart = leaveClusterFrameStart(1:numel(leaveClusterFrameEnd));
            for exitCtr = 1:numel(leaveClusterFrameStart)
                leaveClusterLogInd(leaveClusterFrameStart(exitCtr):leaveClusterFrameEnd(exitCtr))=true;
            end
            leaveClusterLogInd(inClusterLogInd)=false; % exclude when worms move back into a cluster
            loneWormLogInd = min_neighbr_dist>=minNeighbrDist;
            leaveClusterLogInd(loneWormLogInd)=false; % exclude worms that have become lone worm
            leaveClusterLogInd = ismember(features.skeleton_id+1,find(trajData.filtered & leaveClusterLogInd)); % make clusterClusterLogInd the same size as features.filtered
            % find lone worms
            loneWormLogInd = ismember(features.skeleton_id+1,find(trajData.filtered & loneWormLogInd));
            %% calculate or extract desired feature values
            % calculte integrated angular speeds
            uniquePaths = unique(features.worm_index);
            iHAS_leaveCluster_poolPaths = cell(1,numel(uniquePaths));
            iHAS_loneWorm_poolPaths = cell(1,numel(uniquePaths));
            for pathCtr = 2:numel(uniquePaths)
                path = uniquePaths(pathCtr);
                pathLogInd = ismember(features.worm_index, path);
                % leaveCluster
                HAS_leaveClusterLogInd = pathLogInd & features.filtered & leaveClusterLogInd;
                HAS_leaveCluster = features.head_tip_motion_direction(HAS_leaveClusterLogInd);
                if size(HAS_leaveCluster,1)> postExitDuration*frameRate % only concerned about turn events within specified seconds of exiting cluster
                    HAS_leaveCluster = HAS_leaveCluster(1:postExitDuration*frameRate);
                end
                HAS_leaveCluster(isnan(HAS_leaveCluster))=[];
                if ~isempty(HAS_leaveCluster)
                    HAS_leaveClusterContTurns = vertcat(1, find(HAS_leaveCluster(1:end-1).*...
                        HAS_leaveCluster(2:end)<0), numel(HAS_leaveCluster));
                    % generate a list of positions of HAS_leaveCluster vector where signs change
                    numContTurns_leaveCluster = length(HAS_leaveClusterContTurns)-1;
                    iHAS_singlePath = NaN(1,numContTurns_leaveCluster); % create variable to hold integrated values from each continuous turn pattern
                    for contTurnCtr = 1:numContTurns_leaveCluster
                        iHASContTurn = abs(sum(HAS_leaveCluster(HAS_leaveClusterContTurns(contTurnCtr)):...
                            HAS_leaveCluster(HAS_leaveClusterContTurns(contTurnCtr+1)))); % integration
                        iHAS_singlePath(contTurnCtr) = iHASContTurn;
                    end
                    iHAS_leaveCluster_poolPaths{pathCtr}=iHAS_singlePath;
                end
                iHAS_leaveCluster_pathPooled = horzcat(iHAS_leaveCluster_poolPaths{:});
                iHAS_leaveCluster_poolFiles{fileCtr} = iHAS_leaveCluster_pathPooled;
                % lone worm
                HAS_loneWormLogInd = pathLogInd & features.filtered & loneWormLogInd;
                HAS_loneWorm = features.head_tip_motion_direction(HAS_loneWormLogInd);
                HAS_loneWorm(isnan(HAS_loneWorm))=[];
                if ~isempty(HAS_loneWorm) & size(HAS_loneWorm,1)>=postExitDuration*frameRate
                    initialFrame_loneWorm = randi(size(HAS_loneWorm,1)-postExitDuration*frameRate+1,1);
                    HAS_loneWorm = HAS_loneWorm(initialFrame_loneWorm:initialFrame_loneWorm+postExitDuration*frameRate-1); % randomly sample 5s of consecutive event
                    HAS_loneWormContTurns = vertcat(1, find(HAS_loneWorm(1:end-1).*...
                        HAS_loneWorm(2:end)<0), numel(HAS_loneWorm));
                    % generate a list of positions of HAS_loneWorm vector where signs change
                    numContTurns_loneWorm = length(HAS_loneWormContTurns)-1;
                    iHAS_singlePath = NaN(1,numContTurns_loneWorm); % create variable to hold integrated values from each continuous turn pattern
                    for contTurnCtr = 1:numContTurns_loneWorm
                        iHASContTurn = abs(sum(HAS_loneWorm(HAS_loneWormContTurns(contTurnCtr)):...
                            HAS_loneWorm(HAS_loneWormContTurns(contTurnCtr+1)))); % integration
                        iHAS_singlePath(contTurnCtr) = iHASContTurn;
                    end
                    iHAS_loneWorm_poolPaths{pathCtr}=iHAS_singlePath;
                end
                iHAS_loneWorm_pathPooled = horzcat(iHAS_loneWorm_poolPaths{:});
                iHAS_loneWorm_poolFiles{fileCtr} = iHAS_loneWorm_pathPooled;                
            end
        end
        iHAS_leaveCluster_poolFiles = horzcat(iHAS_leaveCluster_poolFiles{:})';
        iHAS_leaveCluster_poolFiles = iHAS_leaveCluster_poolFiles(iHAS_leaveCluster_poolFiles>0);
        iHAS_loneWorm_poolFiles = horzcat(iHAS_loneWorm_poolFiles{:})';
        iHAS_loneWorm_poolFiles = iHAS_loneWorm_poolFiles(iHAS_loneWorm_poolFiles>0);
        set(0,'CurrentFigure',iHASFig)
        histogram(iHAS_leaveCluster_poolFiles, 'Normalization','pdf','DisplayStyle','stairs')
        histogram(iHAS_loneWorm_poolFiles, 'Normalization','pdf','DisplayStyle','stairs')
        title(['intHeadAngularSpeed\_' strains{strainCtr} '\_' wormnums{numCtr}])
        legend('leaveCluster','loneWorms')
        set(iHASFig,'PaperUnits','centimeters')
        figurename = ['figures/turns/intHeadAngularSpeed_' strains{strainCtr} '_' wormnums{numCtr} '_' phase '_CL'];
        exportfig(iHASFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
    end
end