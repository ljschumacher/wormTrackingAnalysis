% plot head angular speed distribution

%% issues to address:
% How can some of the straighter-looking trajectories produce high angles?
% Truncate N2 traj according to npr1 traj length distribution - circumvet issues caused by half turns etc. 
% No leave cluster traj for some npr1 movies
% The traj over 5 seconds look quite short in micron terms compared with worm length

clear
close all

%% set parameters
phase = 'fullMovie'; % 'fullMovie', 'joining', or 'sweeping'.
dataset = 2; % 1 or 2
marker = 'pharynx'; % 'pharynx' or 'bodywall'
strains = {'npr1','N2'}; % {'npr1','N2'}
wormnums = {'40'};% {'40'};
wormcats = {'leaveCluster','loneWorm'}; %'leaveCluster','loneWorm'
smoothing = true;
postExitDuration = 5; % set the duration (in seconds) after a worm exits a cluster to be included in the analysis
headAngSpeedRanges = [0 0.5; 1 1.5; 2 2.5; 3 3.5; 4 20]; % one div
%headAngSpeedRanges = [0 2; 5 7; 10 12; 20 25; 26 50]; % second div
visualiseAngSpeedRangeSamples = true; % true or false
numSampleTraj = 5; % number of sample trajectories to be plotted for each angular speed range
pixelsize = 100/19.5; % 100 microns are 19.5 pixels

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
            [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_list_hamm.xlsx'],1,'A1:E15','basic');
        elseif dataset ==2 & strcmp(marker,'pharynx')
            [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_g_list_hamm.xlsx'],1,'A1:E15','basic');
        elseif dataset ==2 & strcmp(marker,'bodywall')
            [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_r_list_hamm.xlsx'],1,'A1:E15','basic');
        else
            warning('specified dataset/marker combination does not exist')
        end
        numFiles = length(filenames);
        % create empty cell arrays to hold individual file values, so they can be pooled for a given strain/density combination
        for wormcatCtr = 1:length(wormcats)
            headAngSpeed.(wormcats{wormcatCtr}) = cell(numFiles,1);
            frameRunLengths.(wormcats{wormcatCtr}) = cell(numFiles,1);
        end
        
        if visualiseAngSpeedRangeSamples
            % save up to 500 sets of xy coordinates for paths that fall within a certain angular speed range to plot sample trajectories
            for wormcatCtr = 1:length(wormcats)
                headAngSpeedSampleTraj.(wormcats{wormcatCtr}) = cell(500,2,size(headAngSpeedRanges,1));
                headAngSpeedSampleCtr.(wormcats{wormcatCtr}) = ones(size(headAngSpeedRanges,1),1);
            end
        end
        
        %% go through individual movies
        for fileCtr = 1:numFiles
            %% load data
            filename = filenames{fileCtr}
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton'); % in pixels
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            
            %% filter worms by various criteria
            if strcmp(marker, 'pharynx')
                % filter green by blob size and intensity
                trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                    intensityThresholds_g(wormnum),maxBlobSize_g);
            elseif strcmp(marker, 'bodywall')
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
            [leaveClusterLogInd, loneWormLogInd] = findLeaveClusterWorms(filename,inClusterNeighbourNum,minNeighbrDist,postExitDuration);
            
            %% calculate or extract desired feature values
            % obtain all xy coordinates
            worm_xcoords = squeeze(skelData(1,:,:))';
            worm_ycoords = squeeze(skelData(2,:,:))';
            uniqueWorms = unique(trajData.worm_index_joined);
            
            % initialise
            for wormcatCtr = 1:length(wormcats)
                headAngSpeed.(wormcats{wormcatCtr}){fileCtr} = NaN(numel(uniqueWorms),100); % assume each worm has up to 100 leave cluster traj. Checks for this later and warns if not enough
                frameRunLengths.(wormcats{wormcatCtr}){fileCtr} = [];
            end
            
            % loop through each worm path
            for wormCtr = 1:numel(uniqueWorms)
                wormLogInd = trajData.worm_index_joined==uniqueWorms(wormCtr) & trajData.filtered;
                for wormcatCtr = 1:length(wormcats)
                    wormCatLogInd = wormLogInd & eval([wormcats{wormcatCtr} 'LogInd']);
                    wormFrames = trajData.frame_number(wormCatLogInd)';
                    % break down frames into continuous trajectories
                    if ~isempty(wormFrames)
                        continuousFrameRuns = getContinuousTraj(wormFrames);
                        % save the length of each continuous trajectory
                        frameRunLengths.(wormcats{wormcatCtr}){fileCtr} = [frameRunLengths.(wormcats{wormcatCtr}){fileCtr} cellfun(@numel,continuousFrameRuns)];
                        % filter for minimum traj length
                        continuousFramesMinLengthLogInd = cellfun(@numel,continuousFrameRuns)>=frameRate+1;
                        continuousFrameRuns = continuousFrameRuns(continuousFramesMinLengthLogInd);
                        % go through each traj that fits the min length criteria to obtain xy coordinates
                        if size(continuousFrameRuns,2)>100
                            warning('more trajectories present than allocated space to hold values for')
                        end
                        % loop through each trajectory
                        for trajCtr = 1:size(continuousFrameRuns,2)
                            continuousframes = continuousFrameRuns{trajCtr};
                            wormtraj_xcoords = NaN(numel(continuousframes),size(worm_xcoords,2));
                            wormtraj_ycoords = NaN(numel(continuousframes),size(worm_ycoords,2));
                            for trajRunFrameCtr = 1:numel(continuousframes)
                                wormtraj_xcoords(trajRunFrameCtr,:) = worm_xcoords(wormCatLogInd...
                                    & trajData.frame_number == continuousframes(trajRunFrameCtr),:);
                                wormtraj_ycoords(trajRunFrameCtr,:) = worm_ycoords(wormCatLogInd...
                                    & trajData.frame_number == continuousframes(trajRunFrameCtr),:);
                            end
                            % filter for maximum traj length
                            maxTotalFrames = postExitDuration*frameRate;
                            if size(wormtraj_xcoords,1)> maxTotalFrames
                                if strcmp(wormcats{wormcatCtr},'leaveCluster')
                                    % always start leave cluster traj from the moment the worm exits cluster
                                    firstFrame = 1;
                                elseif strcmp(wormcats{wormcatCtr},'loneWorm')
                                    % randomly sample the start of lone worm traj
                                    firstFrame = randi(size(wormtraj_xcoords,1)-maxTotalFrames,1);
                                end
                                % truncate the traj at maximum length
                                lastFrame = firstFrame + maxTotalFrames;
                                wormtraj_xcoords = wormtraj_xcoords(firstFrame:lastFrame,:);
                                wormtraj_ycoords = wormtraj_ycoords(firstFrame:lastFrame,:);
                            end
                            % take mean head coordinates
                            if strcmp(marker,'bodywall')
                                wormhead_xcoords = nanmean(wormtraj_xcoords(:,1:8),2)'; % take mean of the first 8 out of 49 nodes (head only)
                                wormhead_ycoords = nanmean(wormtraj_ycoords(:,1:8),2)';
                            elseif strcmp(marker,'pharynx')
                                wormhead_xcoords = nanmean(wormtraj_xcoords,2)'; % pharynx marker has 2 nodes so take the mean of those 2
                                wormhead_ycoords = nanmean(wormtraj_ycoords,2)';
                            end
                            % calculate head angles
                            [angleArray,meanAngles] = makeAngleArray(wormhead_xcoords,wormhead_ycoords);
                            headAngle = angleArray+meanAngles;
                            % get angle difference
                            if smoothing
                                smoothFactor = frameRate+1; % set smoothing to be over 1 second
                                angleDiff = headAngle(2:end) - headAngle(1:end-1);
                                totalSmoothedFrames = length(angleDiff)-smoothFactor;
                                smoothedHeadAngleDiff = NaN(totalSmoothedFrames,1);
                                for smoothFrameCtr = 1:totalSmoothedFrames
                                    smoothedHeadAngleDiff(smoothFrameCtr) = nanmean(angleDiff(smoothFrameCtr:smoothFrameCtr+smoothFactor));
                                end
                                % calculate total smoothed head angle change per second
                                headAngSpeed.(wormcats{wormcatCtr}){fileCtr}(wormCtr,trajCtr) =...
                                    abs(nansum(smoothedHeadAngleDiff)/totalSmoothedFrames*frameRate);
                            else
                                headAngleDiff = headAngle(2:end) - headAngle(1:end-1);
                                % calculate total smoothed head angle change per second
                                headAngSpeed.(wormcats{wormcatCtr}){fileCtr}(wormCtr,trajCtr) =...
                                    abs(nansum(headAngleDiff)/length(headAngle)*frameRate);
                            end
                            % save xy coordinates for paths that fall within certain angular speed ranges
                            if visualiseAngSpeedRangeSamples
                                % loop through each range to see which one it falls within
                                for rangeCtr = 1:size(headAngSpeedRanges,1)
                                    if headAngSpeedRanges(rangeCtr,1)<headAngSpeed.(wormcats{wormcatCtr}){fileCtr}(wormCtr,trajCtr) & ...
                                        headAngSpeed.(wormcats{wormcatCtr}){fileCtr}(wormCtr,trajCtr)<headAngSpeedRanges(rangeCtr,2)
                                        % save xy coordinates in microns
                                        wormhead_xcoords = wormhead_xcoords * pixelsize; % turn pixels into microns
                                        wormhead_ycoords = wormhead_ycoords * pixelsize;
                                        headAngSpeedSampleTraj.(wormcats{wormcatCtr}){headAngSpeedSampleCtr.(wormcats{wormcatCtr})(rangeCtr),1,rangeCtr} = wormhead_xcoords;
                                        headAngSpeedSampleTraj.(wormcats{wormcatCtr}){headAngSpeedSampleCtr.(wormcats{wormcatCtr})(rangeCtr),2,rangeCtr} = wormhead_ycoords;
                                        % update the traj counter for that range (until up to 500 preallocated trajectory spaces)
                                        if headAngSpeedSampleCtr.(wormcats{wormcatCtr})(rangeCtr) < size(headAngSpeedSampleTraj.(wormcats{wormcatCtr}),1)
                                           headAngSpeedSampleCtr.(wormcats{wormcatCtr})(rangeCtr) = headAngSpeedSampleCtr.(wormcats{wormcatCtr})(rangeCtr)+1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            for wormcatCtr = 1:length(wormcats)
                headAngSpeed.(wormcats{wormcatCtr}){fileCtr}(isnan(headAngSpeed.(wormcats{wormcatCtr}){fileCtr}))=[];
            end
        end
        % pool data from all files belonging to the same strain and worm density
        for wormcatCtr = 1:length(wormcats)
            headAngSpeedPool.(wormcats{wormcatCtr}) = horzcat(headAngSpeed.(wormcats{wormcatCtr}){:});
            frameRunLengths.(wormcats{wormcatCtr}) = horzcat(frameRunLengths.(wormcats{wormcatCtr}){:});
        end
        
        %% plot data, format, and export
        % save angular speed values
        %save(['figures/turns/results/headAngSpeed_' strains{strainCtr} '_' wormnums{numCtr} '.mat'],'headAngSpeed')
        % plot angular speed distribution
        headAngSpeedFig = figure; hold on
        for wormcatCtr = 1:length(wormcats)
            histogram(headAngSpeedPool.(wormcats{wormcatCtr}),'Normalization','pdf','DisplayStyle','stairs')
            Legend{wormcatCtr} = strcat(wormcats{wormcatCtr}, ',n=', num2str(size(headAngSpeedPool.(wormcats{wormcatCtr}),2)));
        end
        legend(Legend)
        title([strains{strainCtr} '\_' wormnum],'FontWeight','normal')
        xlabel('head angle change rate (°/s)')
        ylabel('probability')
        xlim([0 12])
        ylim([0 2])
        set(headAngSpeedFig,'PaperUnits','centimeters')
        figurename = ['figures/turns/headAngSpeed_' strain '_' wormnum '_' phase '_data' num2str(dataset) '_' marker];
        %savefig(headAngSpeedFig,[figurename '.fig'])
        load('exportOptions.mat')
        %exportfig(headAngSpeedFig,[figurename '.eps'],exportOptions)
        %system(['epstopdf ' figurename '.eps']);
        %system(['rm ' figurename '.eps']);
        
        if visualiseAngSpeedRangeSamples
            %save trajectories
            save(['figures/turns/results/headAngSpeedSampleTraj_' strain '_' wormnum '.mat'],'headAngSpeedSampleTraj')
            %plot trajectories
            plotSampleHeadAngSpeedTraj(headAngSpeedRanges,wormcats,numSampleTraj,strain,wormnum);
        end
    end
end