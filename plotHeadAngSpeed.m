% plot cumulative smoothed head angle changes

clear
close all

%% set parameters
phase = 'fullMovie'; % 'fullMovie', 'joining', or 'sweeping'.
dataset = 2; % 1 or 2
marker = 'pharynx'; % 'pharynx' or 'bodywall'
strains = {'npr1','N2'}; % {'npr1','N2'}
wormnums = {'40'};% {'40'};
wormcat = {'leaveCluster','loneWorm'}; %'leaveCluster','loneWorm'
postExitDuration = 5; % set the duration (in seconds) after a worm exits a cluster to be included in the analysis
headAngSpeedRanges = [0 2; 10 12; 20 22; 30 32; 50 52];
%headAngSpeedRanges = [0 1; 6 8; 13 18; 21 23; 27 30]; % ranges in degrees/s for plotting sample trajectories
visualiseAngSpeedRangeSamples = true; % true or false

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
        for wormcatCtr = 1:length(wormcat)
            headAngSpeed.(wormcat{wormcatCtr}) = cell(numFiles,1);
            frameRunLengths.(wormcat{wormcatCtr}) = cell(numFiles,1);
        end
        
        if visualiseAngSpeedRangeSamples
            % save up to 500 sets of xy coordinates for paths that fall within a certain angular speed range to plot sample trajectories
            for wormcatCtr = 1:length(wormcat)
                headAngleSpeedSampleTraj.(wormcat{wormcatCtr}) = cell(500,2,size(headAngSpeedRanges,1));
                headAngleSpeedSampleCtr.(wormcat{wormcatCtr}) = ones(size(headAngSpeedRanges,1),1);
            end
        end
        
        %% go through individual movies
        for fileCtr = 1:numFiles
            %% load data
            filename = filenames{fileCtr}
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton');
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            smoothFactor = 1*frameRate; % set smoothing to be over 1 second
            
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
            % restrict to only forward-moving worms
            signedSpeedLogInd = blobFeats.signed_speed>=0;
            trajData.filtered(~signedSpeedLogInd) = false;
            % find worms that have just left a cluster vs lone worms
            [leaveClusterLogInd, loneWormLogInd] = findLeaveClusterWorms(filename,inClusterNeighbourNum,minNeighbrDist,postExitDuration);
            
            %% calculate or extract desired feature values
            % obtain all xy coordinates
            worm_xcoords = squeeze(skelData(1,:,:))';
            worm_ycoords = squeeze(skelData(2,:,:))';
            uniqueWorms = unique(trajData.worm_index_joined);
            
            % initialise
            for wormcatCtr = 1:length(wormcat)
                headAngSpeed.(wormcat{wormcatCtr}){fileCtr} = NaN(numel(uniqueWorms),100); % assume each worm has up to 100 leave cluster traj. Checks for this later and warns if not enough
                frameRunLengths.(wormcat{wormcatCtr}){fileCtr} = [];
            end
            
            % loop through each worm path
            for wormCtr = 1:numel(uniqueWorms)
                wormLogInd = trajData.worm_index_joined==uniqueWorms(wormCtr) & trajData.filtered;
                for wormcatCtr = 1:length(wormcat)
                    wormCatLogInd = wormLogInd & eval([wormcat{wormcatCtr} 'LogInd']);
                    wormFrames = trajData.frame_number(wormCatLogInd)';
                    % break down frames into continuous trajectories
                    if ~isempty(wormFrames)
                        continuousFrameRuns = getContinuousTraj(wormFrames);
                        % save the length of each continuous trajectory
                        frameRunLengths.(wormcat{wormcatCtr}){fileCtr} = [frameRunLengths.(wormcat{wormcatCtr}){fileCtr} cellfun(@numel,continuousFrameRuns)];
                        % filter for minimum traj length
                        continuousFramesMinLengthLogInd = cellfun(@numel,continuousFrameRuns)>=smoothFactor;
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
                            if size(wormtraj_xcoords,1)>(postExitDuration)*frameRate
                                if strcmp(wormcat{wormcatCtr},'leaveCluster')
                                    % always start leave cluster traj from the moment the worm exits cluster
                                    firstFrame = 1;
                                elseif strcmp(wormcat{wormcatCtr},'loneWorm')
                                    % randomly sample the start of lone worm traj
                                    firstFrame = randi(size(wormtraj_xcoords,1)-postExitDuration*frameRate,1);
                                end
                                % truncate the traj at maximum length
                                lastFrame = firstFrame + postExitDuration*frameRate;
                                wormtraj_xcoords = wormtraj_xcoords(firstFrame:lastFrame,:);
                                wormtraj_ycoords = wormtraj_ycoords(firstFrame:lastFrame,:);
                            end
                            % calculate angles
                            [angleArray,meanAngles] = makeAngleArray(wormtraj_xcoords,wormtraj_ycoords);
                            angleArray = angleArray+meanAngles;
                            % take mean head angles
                            if strcmp(marker,'bodywall')
                                headAngle = nanmean(angleArray(:,1:8),2); % take mean of the first 8 out of 49 nodes (head only)
                            elseif strcmp(marker,'pharynx')
                                headAngle = nanmean(angleArray,2); % pharynx marker has 2 nodes so take the mean of those 2
                            end
                            % get angle difference over 1 second (over 1 second rather than between each frame to implement smoothing)
                            headAngleDiff = headAngle(smoothFactor:end) - headAngle(1:end-smoothFactor+1);
                            % calculate total smoothed head angle change per second
                            headAngSpeed.(wormcat{wormcatCtr}){fileCtr}(wormCtr,trajCtr) =...
                                abs(nansum(headAngleDiff)/length(headAngleDiff)*frameRate);
                            % save xy coordinates for paths that fall within certain angular speed ranges
                            if visualiseAngSpeedRangeSamples
                                % loop through each range to see which one it falls within
                                for rangeCtr = 1:size(headAngSpeedRanges,1)
                                    if headAngSpeedRanges(rangeCtr,1)<headAngSpeed.(wormcat{wormcatCtr}){fileCtr}(wormCtr,trajCtr) & ...
                                            headAngSpeed.(wormcat{wormcatCtr}){fileCtr}(wormCtr,trajCtr)<headAngSpeedRanges(rangeCtr,2)
                                        % save xy coordinates
                                        headAngleSpeedSampleTraj.(wormcat{wormcatCtr}){headAngleSpeedSampleCtr.(wormcat{wormcatCtr})(rangeCtr),1,rangeCtr} = wormtraj_xcoords;
                                        headAngleSpeedSampleTraj.(wormcat{wormcatCtr}){headAngleSpeedSampleCtr.(wormcat{wormcatCtr})(rangeCtr),2,rangeCtr} = wormtraj_ycoords;
                                        % update the traj counter for that range (until up to 500 preallocated trajectory spaces)
                                        if headAngleSpeedSampleCtr.(wormcat{wormcatCtr})(rangeCtr) < size(headAngleSpeedSampleTraj.(wormcat{wormcatCtr}),1)
                                            headAngleSpeedSampleCtr.(wormcat{wormcatCtr})(rangeCtr) = headAngleSpeedSampleCtr.(wormcat{wormcatCtr})(rangeCtr)+1;
                                        end
                                        % remove empty rows if less than 500 traj are saved
                                        if headAngleSpeedSampleCtr.(wormcat{wormcatCtr})(rangeCtr)<500
                                           headAngleSpeedSampleTraj.(wormcat{wormcatCtr}) = headAngleSpeedSampleTraj.(wormcat{wormcatCtr})(headAngleSpeedSampleCtr.(wormcat{wormcatCtr})(rangeCtr),2,(rangeCtr));
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            for wormcatCtr = 1:length(wormcat)
                headAngSpeed.(wormcat{wormcatCtr}){fileCtr}(isnan(headAngSpeed.(wormcat{wormcatCtr}){fileCtr}))=[];
            end
        end
        % pool data from all files belonging to the same strain and worm density
        for wormcatCtr = 1:length(wormcat)
            headAngSpeedPool.(wormcat{wormcatCtr}) = horzcat(headAngSpeed.(wormcat{wormcatCtr}){:});
            frameRunLengths.(wormcat{wormcatCtr}) = horzcat(frameRunLengths.(wormcat{wormcatCtr}){:});
        end
        
        %% plot data, format, and export
        % save angular speed values
        save(['figures/turns/results/headAngSpeed_' strains{strainCtr} '_' wormnums{numCtr} '.mat'],'headAngSpeed')
        % plot angular speed distribution
        headAngSpeedFig = figure; hold on
        for wormcatCtr = 1:length(wormcat)
            histogram(headAngSpeedPool.(wormcat{wormcatCtr}),'Normalization','pdf','DisplayStyle','stairs')
            Legend{wormcatCtr} = strcat(wormcat{wormcatCtr}, ',n=', num2str(size(headAngSpeed.(wormcat{wormcatCtr}),2)));
        end
        legend(Legend)
        title([strains{strainCtr} '\_' wormnums{numCtr}],'FontWeight','normal')
        xlabel('head angle change rate (°/s)')
        ylabel('probability')
        %xlim([0 30])
        %ylim([0 0.09])
        set(headAngSpeedFig,'PaperUnits','centimeters')
        figurename = ['figures/turns/headAngSpeed_' strains{strainCtr} '_' wormnums{numCtr} '_' phase '_data' num2str(dataset) '_' marker];
        savefig(headAngSpeedFig,[figurename '.fig'])
        load('exportOptions.mat')
        exportfig(headAngSpeedFig,[figurename '.eps'],exportOptions)
        %system(['epstopdf ' figurename '.eps']);
        %system(['rm ' figurename '.eps']);
        
        if visualiseAngSpeedRangeSamples
            %save trajectories
            save(['figures/turns/results/headAngSpeedSampleTraj_' strains{strainCtr} '_' wormnums{numCtr} '.mat'],'headAngleSpeedSampleTraj')
            %plot trajectories
            numSampleTraj = 5;
            plotSampleHeadAngSpeedTraj(headAngSpeedRanges,wormcat,strain,numSampleTraj);
        end
    end
end