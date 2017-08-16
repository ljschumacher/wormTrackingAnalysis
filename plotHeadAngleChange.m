% plot cumulative smoothed head angle changes

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
phase = 'fullMovie'; % 'fullMovie', 'joining', or 'sweeping'.
dataset = 2; % 1 or 2
marker = 'pharynx'; % 'pharynx' or 'bodywall'
strains = {'npr1','N2'}; % {'npr1','N2'}
wormnums = {'40'};% {'40'};
postExitDuration = 5; % set the duration (in seconds) after a worm exits a cluster to be included in the analysis
sHeadAngSpeedRanges = [0 1; 6 8; 13 18; 21 23; 27 30]; % ranges in degrees/s for plotting sample trajectories

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
        phaseFrames = phaseFrames-1; % correct for python indexing at 0
        numFiles = length(filenames);
        % create empty cell arrays to hold individual file values, so they can be pooled for a given strain/density combination
        sHeadAngleChangeRate_leaveCluster = cell(numFiles,1);
        sHeadAngleChangeRate_loneWorm = cell(numFiles,1);
        headAngleChangeRate_leaveCluster = cell(numFiles,1);
        headAngleChangeRate_loneWorm = cell(numFiles,1);
        sHeadAngleChangeRateFig = figure;
        
        % save up to 500 sets of xy coordinates for paths that fall within a certain angular speed range to plot sample trajectories
        headAngSpeedSample_leaveCluster1 = cell(500,2);
        headAngSpeedSample_leaveCluster2 = cell(500,2);
        headAngSpeedSample_leaveCluster3 = cell(500,2);
        headAngSpeedSample_leaveCluster4 = cell(500,2);
        headAngSpeedSample_leaveCluster5 = cell(500,2);
        sHeadAngleSpeedSample_leaveCluster1 = cell(500,1);
        sHeadAngleSpeedSample_leaveCluster2 = cell(500,1);
        sHeadAngleSpeedSample_leaveCluster3 = cell(500,1);
        sHeadAngleSpeedSample_leaveCluster4 = cell(500,1);
        sHeadAngleSpeedSample_leaveCluster5 = cell(500,1);
        headAngSpeedSample_loneWorm1 = cell(500,2);
        headAngSpeedSample_loneWorm2 = cell(500,2);
        headAngSpeedSample_loneWorm3 = cell(500,2);
        headAngSpeedSample_loneWorm4 = cell(500,2);
        headAngSpeedSample_loneWorm5 = cell(500,2);
        sHeadAngleSpeedSample_loneWorm1 = cell(500,1);
        sHeadAngleSpeedSample_loneWorm2 = cell(500,1);
        sHeadAngleSpeedSample_loneWorm3 = cell(500,1);
        sHeadAngleSpeedSample_loneWorm4 = cell(500,1);
        sHeadAngleSpeedSample_loneWorm5 = cell(500,1);
        headAngSpeedSample_leaveCluster1Ctr = 1;
        headAngSpeedSample_leaveCluster2Ctr = 1;
        headAngSpeedSample_leaveCluster3Ctr = 1;
        headAngSpeedSample_leaveCluster4Ctr = 1;
        headAngSpeedSample_leaveCluster5Ctr = 1;
        headAngSpeedSample_loneWorm1Ctr = 1;
        headAngSpeedSample_loneWorm2Ctr = 1;
        headAngSpeedSample_loneWorm3Ctr = 1;
        headAngSpeedSample_loneWorm4Ctr = 1;
        headAngSpeedSample_loneWorm5Ctr = 1;
        HeadAngSpeedSample_leaveCluster1Fig = figure; hold on
        HeadAngSpeedSample_leaveCluster2Fig = figure; hold on
        HeadAngSpeedSample_leaveCluster3Fig = figure; hold on
        HeadAngSpeedSample_leaveCluster4Fig = figure; hold on
        HeadAngSpeedSample_leaveCluster5Fig = figure; hold on
        HeadAngSpeedSample_loneWorm1Fig = figure; hold on
        HeadAngSpeedSample_loneWorm2Fig = figure; hold on
        HeadAngSpeedSample_loneWorm3Fig = figure; hold on
        HeadAngSpeedSample_loneWorm4Fig = figure; hold on
        HeadAngSpeedSample_loneWorm5Fig = figure; hold on
        frameRunLengths_leaveCluster = cell(numFiles,1);
        frameRunLengths_loneWorm = cell(numFiles,1);
        
        %% go through individual movies
        for fileCtr = 1:numFiles
            %% load data
            filename = filenames{fileCtr}
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton');
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            smoothFactor = 1*frameRate; % smooth head angles over 1 second
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
            phaseFrameLogInd = trajData.frame_number <= lastFrame & trajData.frame_number >= firstFrame;
            trajData.filtered(~phaseFrameLogInd) = false;
            % find worms that have just left a cluster vs lone worms
            [leaveClusterLogInd, loneWormLogInd] = findLeaveClusterWorms(filename,inClusterNeighbourNum,minNeighbrDist,postExitDuration);
            
            %% calculate or extract desired feature values
            % calculate cHeadAngle by looping through each worm path
            uniqueWorms = unique(trajData.worm_index_joined);
            worm_xcoords = squeeze(skelData(1,:,:))';
            worm_ycoords = squeeze(skelData(2,:,:))';
            % initialise
            headAngleChangeRate_leaveCluster_thisFile = NaN(numel(uniqueWorms),100);
            headAngleChangeRate_loneWorm_thisFile = NaN(numel(uniqueWorms),100);
            sHeadAngleChangeRate_leaveCluster_thisFile = NaN(numel(uniqueWorms),100); % assume each worm has up to 50 leave cluster traj. Checks for this later and warns if not enough
            sHeadAngleChangeRate_loneWorm_thisFile = NaN(numel(uniqueWorms),1000);
            frameRunLengths_leaveCluster_thisFile = [];
            frameRunLengths_loneWorm_thisFile = [];
            
            for wormpathCtr = 1:numel(uniqueWorms)
                wormpathLogInd = trajData.worm_index_joined==uniqueWorms(wormpathCtr) & trajData.filtered;
                wormpathFrames_leaveCluster = trajData.frame_number(wormpathLogInd & leaveClusterLogInd)';
                wormpathFrames_loneWorm = trajData.frame_number(wormpathLogInd & loneWormLogInd)';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % break down leave cluster frames into continuous trajectories
                if ~isempty(wormpathFrames_leaveCluster)
                    frameStepSize = diff(wormpathFrames_leaveCluster);
                    continuousFrameRuns = {};
                    continuousRunCtr = 1;
                    continuousFrameRuns{continuousRunCtr} = wormpathFrames_leaveCluster(1);
                    for i = 1:numel(frameStepSize)
                        if frameStepSize(i) == 1
                            continuousFrameRuns{continuousRunCtr} = [continuousFrameRuns{continuousRunCtr} wormpathFrames_leaveCluster(i+1)];
                        else
                            continuousRunCtr = continuousRunCtr + 1;
                            continuousFrameRuns{continuousRunCtr} = wormpathFrames_leaveCluster(i+1);
                        end
                    end
                    % filter for minimum traj length
                    frameRunLengths_leaveCluster_thisFile =  [frameRunLengths_leaveCluster_thisFile cellfun(@numel,continuousFrameRuns)];
                    continuousFramesMinLengthLogInd = cellfun(@numel,continuousFrameRuns)>=round(1.5*frameRate);
                    continuousFrameRuns = continuousFrameRuns(continuousFramesMinLengthLogInd); % gives cell arrays containing frame numbers for the chosen worm where frames are continuous
                    % go through each traj that fits the min length criteria to obtain xy coordinates
                    if size(continuousFrameRuns,2)>100
                        warning('more trajectories present than allocated space to hold values for')
                    end
                    for trajRunCtr = 1:size(continuousFrameRuns,2)
                        continuousframes = continuousFrameRuns{trajRunCtr};
                        wormpath_xcoords_leaveCluster = NaN(numel(continuousframes),size(worm_xcoords,2));
                        wormpath_ycoords_leaveCluster = NaN(numel(continuousframes),size(worm_ycoords,2));
                        for trajRunFrameCtr = 1:numel(continuousframes)
                            wormpath_xcoords_leaveCluster(trajRunFrameCtr,:) = worm_xcoords((wormpathLogInd & leaveClusterLogInd...
                                & trajData.frame_number == continuousframes(trajRunFrameCtr)),:);
                            wormpath_ycoords_leaveCluster(trajRunFrameCtr,:) = worm_ycoords((wormpathLogInd & leaveClusterLogInd...
                                & trajData.frame_number == continuousframes(trajRunFrameCtr)),:);
                        end
                        % filter for maximum traj length
                        if size(wormpath_xcoords_leaveCluster,1)>(postExitDuration+1)*frameRate
                            firstFrame = 1;
                            lastFrame = firstFrame + (postExitDuration+1)*frameRate;
                            wormpath_xcoords_leaveCluster = wormpath_xcoords_leaveCluster(firstFrame:lastFrame,:);
                            wormpath_ycoords_leaveCluster = wormpath_ycoords_leaveCluster(firstFrame:lastFrame,:);
                        end
                        % calculate angles
                        [angleArray_leaveCluster,meanAngles_leaveCluster] = makeAngleArray(wormpath_xcoords_leaveCluster,wormpath_ycoords_leaveCluster);
                        angleArray_leaveCluster = angleArray_leaveCluster+meanAngles_leaveCluster;
                        % take mean head angles
                        if strcmp(marker,'bodywall')
                            headAngle_leaveCluster = nanmean(angleArray_leaveCluster(:,1:8),2);
                        elseif strcmp(marker,'pharynx')
                            headAngle_leaveCluster = angleArray_leaveCluster;
                        end
                        % set to head angles to smooth over 1 second
                        smoothHeadAngle_leaveCluster = NaN(length(headAngle_leaveCluster)-smoothFactor,1);
                        for smoothCtr = 1:(length(headAngle_leaveCluster)-smoothFactor)
                            smoothHeadAngle_leaveCluster(smoothCtr) = nanmean(headAngle_leaveCluster(smoothCtr:smoothCtr+smoothFactor));
                        end
                        % calculate total smoothed head angle change per second
                        headAngleChangeRate_leaveCluster_thisFile(wormpathCtr,trajRunCtr) = abs(nansum(headAngle_leaveCluster)/length(headAngle_leaveCluster)/length(headAngle_leaveCluster)*frameRate*frameRate);
                        sHeadAngleChangeRate_leaveCluster_thisFile(wormpathCtr,trajRunCtr) = abs(nansum(smoothHeadAngle_leaveCluster)/length(smoothHeadAngle_leaveCluster)/length(smoothHeadAngle_leaveCluster)*frameRate*frameRate);
                        % save/plot trajectories that fall within certain angular speed ranges 
                        wormpath_xcoords_leaveCluster = mean(wormpath_xcoords_leaveCluster,2);
                        wormpath_ycoords_leaveCluster = mean(wormpath_ycoords_leaveCluster,2);
                        if sHeadAngSpeedRanges(1,1)<sHeadAngleChangeRate_leaveCluster_thisFile(wormpathCtr,trajRunCtr) & ...
                                sHeadAngleChangeRate_leaveCluster_thisFile(wormpathCtr,trajRunCtr)<sHeadAngSpeedRanges(1,2)
                            headAngSpeedSample_leaveCluster1{headAngSpeedSample_leaveCluster1Ctr,1} = wormpath_xcoords_leaveCluster;
                            headAngSpeedSample_leaveCluster1{headAngSpeedSample_leaveCluster1Ctr,2} = wormpath_ycoords_leaveCluster;
                            sHeadAngleSpeedSample_leaveCluster1{headAngSpeedSample_leaveCluster1Ctr} = smoothHeadAngle_leaveCluster;
                            if headAngSpeedSample_leaveCluster1Ctr < size(headAngSpeedSample_leaveCluster1,1)
                                headAngSpeedSample_leaveCluster1Ctr = headAngSpeedSample_leaveCluster1Ctr+1;
                            end
                            set(0,'CurrentFigure',HeadAngSpeedSample_leaveCluster1Fig)
                            plot(wormpath_xcoords_leaveCluster,wormpath_ycoords_leaveCluster)
                        elseif sHeadAngSpeedRanges(2,1)<sHeadAngleChangeRate_leaveCluster_thisFile(wormpathCtr,trajRunCtr) &...
                                sHeadAngleChangeRate_leaveCluster_thisFile(wormpathCtr,trajRunCtr)<sHeadAngSpeedRanges(2,2)
                            headAngSpeedSample_leaveCluster2{headAngSpeedSample_leaveCluster2Ctr,1} = wormpath_xcoords_leaveCluster;
                            headAngSpeedSample_leaveCluster2{headAngSpeedSample_leaveCluster2Ctr,2} = wormpath_ycoords_leaveCluster;
                            sHeadAngleSpeedSample_leaveCluster2{headAngSpeedSample_leaveCluster2Ctr} = smoothHeadAngle_leaveCluster;
                            if headAngSpeedSample_leaveCluster2Ctr < size(headAngSpeedSample_leaveCluster2,1)
                                headAngSpeedSample_leaveCluster2Ctr = headAngSpeedSample_leaveCluster2Ctr+1;
                            end
                            set(0,'CurrentFigure',HeadAngSpeedSample_leaveCluster2Fig)
                            plot(wormpath_xcoords_leaveCluster,wormpath_ycoords_leaveCluster)
                        elseif sHeadAngSpeedRanges(3,1)<sHeadAngleChangeRate_leaveCluster_thisFile(wormpathCtr,trajRunCtr) &...
                                sHeadAngleChangeRate_leaveCluster_thisFile(wormpathCtr,trajRunCtr)<sHeadAngSpeedRanges(3,2)
                            headAngSpeedSample_leaveCluster3{headAngSpeedSample_leaveCluster3Ctr,1} = wormpath_xcoords_leaveCluster;
                            headAngSpeedSample_leaveCluster3{headAngSpeedSample_leaveCluster3Ctr,2} = wormpath_ycoords_leaveCluster;
                            sHeadAngleSpeedSample_leaveCluster3{headAngSpeedSample_leaveCluster3Ctr} = smoothHeadAngle_leaveCluster;
                            if headAngSpeedSample_leaveCluster3Ctr < size(headAngSpeedSample_leaveCluster3,1)
                                headAngSpeedSample_leaveCluster3Ctr = headAngSpeedSample_leaveCluster3Ctr+1;
                            end
                            set(0,'CurrentFigure',HeadAngSpeedSample_leaveCluster3Fig)
                            plot(wormpath_xcoords_leaveCluster,wormpath_ycoords_leaveCluster)
                        elseif sHeadAngSpeedRanges(4,1)<sHeadAngleChangeRate_leaveCluster_thisFile(wormpathCtr,trajRunCtr) &...
                                sHeadAngleChangeRate_leaveCluster_thisFile(wormpathCtr,trajRunCtr)<sHeadAngSpeedRanges(4,2)
                            headAngSpeedSample_leaveCluster4{headAngSpeedSample_leaveCluster4Ctr,1} = wormpath_xcoords_leaveCluster;
                            headAngSpeedSample_leaveCluster4{headAngSpeedSample_leaveCluster4Ctr,2} = wormpath_ycoords_leaveCluster;
                            sHeadAngleSpeedSample_leaveCluster4{headAngSpeedSample_leaveCluster4Ctr} = smoothHeadAngle_leaveCluster;
                            if headAngSpeedSample_leaveCluster4Ctr < size(headAngSpeedSample_leaveCluster4,1)
                                headAngSpeedSample_leaveCluster4Ctr = headAngSpeedSample_leaveCluster4Ctr+1;
                            end
                            set(0,'CurrentFigure',HeadAngSpeedSample_leaveCluster4Fig)
                            plot(wormpath_xcoords_leaveCluster,wormpath_ycoords_leaveCluster)
                        elseif sHeadAngSpeedRanges(5,1)<sHeadAngleChangeRate_leaveCluster_thisFile(wormpathCtr,trajRunCtr) &...
                                sHeadAngleChangeRate_leaveCluster_thisFile(wormpathCtr,trajRunCtr)<sHeadAngSpeedRanges(5,2)
                            headAngSpeedSample_leaveCluster5{headAngSpeedSample_leaveCluster5Ctr,1} = wormpath_xcoords_leaveCluster;
                            headAngSpeedSample_leaveCluster5{headAngSpeedSample_leaveCluster5Ctr,2} = wormpath_ycoords_leaveCluster;
                            sHeadAngleSpeedSample_leaveCluster5{headAngSpeedSample_leaveCluster5Ctr} = smoothHeadAngle_leaveCluster;
                            if headAngSpeedSample_leaveCluster5Ctr < size(headAngSpeedSample_leaveCluster5,1)
                                headAngSpeedSample_leaveCluster5Ctr = headAngSpeedSample_leaveCluster5Ctr+1;
                            end
                            set(0,'CurrentFigure',HeadAngSpeedSample_leaveCluster5Fig)
                            plot(wormpath_xcoords_leaveCluster,wormpath_ycoords_leaveCluster)
                        end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % break down lone worm frames into continuous trajectories
                if ~isempty(wormpathFrames_loneWorm)
                    frameStepSize = diff(wormpathFrames_loneWorm);
                    continuousFrameRuns = {};
                    continuousRunCtr = 1;
                    continuousFrameRuns{continuousRunCtr} = wormpathFrames_loneWorm(1);
                    for i = 1:numel(frameStepSize)
                        if frameStepSize(i) == 1
                            continuousFrameRuns{continuousRunCtr} = [continuousFrameRuns{continuousRunCtr} wormpathFrames_loneWorm(i+1)];
                        else
                            continuousRunCtr = continuousRunCtr + 1;
                            continuousFrameRuns{continuousRunCtr} = wormpathFrames_loneWorm(i+1);
                        end
                    end
                    % filter for minimum traj length
                    if size(continuousFrameRuns,2)>1000
                        warning('more trajectories present than allocated space to hold values for')
                    end
                    frameRunLengths_loneWorm_thisFile =  [frameRunLengths_loneWorm_thisFile cellfun(@numel,continuousFrameRuns)];
                    continuousFramesMinLengthLogInd = cellfun(@numel,continuousFrameRuns)>=round(1.5*frameRate);
                    continuousFrameRuns = continuousFrameRuns(continuousFramesMinLengthLogInd); % gives cell arrays containing frame numbers for the chosen worm where frames are continuous
                    for trajRunCtr = 1:size(continuousFrameRuns,2)
                        continuousframes = continuousFrameRuns{trajRunCtr};
                        wormpath_xcoords_loneWorm = NaN(numel(continuousframes),size(worm_xcoords,2));
                        wormpath_ycoords_loneWorm = NaN(numel(continuousframes),size(worm_ycoords,2));
                        for trajRunFrameCtr = 1:numel(continuousframes)
                            wormpath_xcoords_loneWorm(trajRunFrameCtr,:) = worm_xcoords((wormpathLogInd & loneWormLogInd...
                                & trajData.frame_number == continuousframes(trajRunFrameCtr)),:);
                            wormpath_ycoords_loneWorm(trajRunFrameCtr,:) = worm_ycoords((wormpathLogInd & loneWormLogInd...
                                & trajData.frame_number == continuousframes(trajRunFrameCtr)),:);
                        end
                        %randomly sample the start of lone worm traj if exceeding specified maximum traj duration
                        if size(wormpath_xcoords_loneWorm,1)>(postExitDuration+1)*frameRate
                            firstFrame = randi(size(wormpath_xcoords_loneWorm,1)-(postExitDuration+1)*frameRate,1);
                            lastFrame = firstFrame + (postExitDuration+1)*frameRate;
                            wormpath_xcoords_loneWorm = wormpath_xcoords_loneWorm(firstFrame:lastFrame,:);
                            wormpath_ycoords_loneWorm = wormpath_ycoords_loneWorm(firstFrame:lastFrame,:);
                        end
                        % calculate angles
                        [angleArray_loneWorm,meanAngles_loneWorm] = makeAngleArray(wormpath_xcoords_loneWorm,wormpath_ycoords_loneWorm);
                        angleArray_loneWorm = angleArray_loneWorm + meanAngles_loneWorm;
                        % take mean head angles
                        if strcmp(marker,'bodywall')
                            headAngle_loneWorm = nanmean(angleArray_loneWorm(:,1:8),2);
                        elseif strcmp(marker,'pharynx')
                            headAngle_loneWorm = angleArray_loneWorm;
                        end
                        % set to head angles to smooth over 1 second
                        smoothHeadAngle_loneWorm = NaN(length(headAngle_loneWorm)-smoothFactor,1);
                        for smoothCtr = 1:(length(headAngle_loneWorm)-smoothFactor)
                            smoothHeadAngle_loneWorm(smoothCtr) = nanmean(headAngle_loneWorm(smoothCtr:smoothCtr+smoothFactor));
                        end
                        % calculate total smoothed head angle change per second
                        headAngleChangeRate_loneWorm_thisFile(wormpathCtr,trajRunCtr) = abs(nansum(headAngle_loneWorm)/length(headAngle_loneWorm)/length(headAngle_loneWorm)*frameRate);
                        sHeadAngleChangeRate_loneWorm_thisFile(wormpathCtr,trajRunCtr) = abs(nansum(smoothHeadAngle_loneWorm)/length(smoothHeadAngle_loneWorm)/length(headAngle_loneWorm)*frameRate);
                        % save/plot trajectories that fit certain angular speed range criteria
                        wormpath_xcoords_loneWorm = mean(wormpath_xcoords_loneWorm,2);
                        wormpath_ycoords_loneWorm = mean(wormpath_ycoords_loneWorm,2);
                        if sHeadAngSpeedRanges(1,1)<sHeadAngleChangeRate_loneWorm_thisFile(wormpathCtr,trajRunCtr) &...
                                sHeadAngleChangeRate_loneWorm_thisFile(wormpathCtr,trajRunCtr)<sHeadAngSpeedRanges(1,2)
                            headAngSpeedSample_loneWorm1{headAngSpeedSample_loneWorm1Ctr,1} = wormpath_xcoords_loneWorm;
                            headAngSpeedSample_loneWorm1{headAngSpeedSample_loneWorm1Ctr,2} = wormpath_ycoords_loneWorm;
                            sHeadAngleSpeedSample_loneWorm1{headAngSpeedSample_loneWorm1Ctr} = smoothHeadAngle_loneWorm;
                            if headAngSpeedSample_loneWorm1Ctr < size(headAngSpeedSample_loneWorm1,1)
                                headAngSpeedSample_loneWorm1Ctr = headAngSpeedSample_loneWorm1Ctr+1;
                            end
                            set(0,'CurrentFigure',HeadAngSpeedSample_loneWorm1Fig)
                            plot(wormpath_xcoords_loneWorm,wormpath_ycoords_loneWorm)
                        elseif sHeadAngSpeedRanges(2,1)<sHeadAngleChangeRate_loneWorm_thisFile(wormpathCtr,trajRunCtr) &...
                                sHeadAngleChangeRate_loneWorm_thisFile(wormpathCtr,trajRunCtr)<sHeadAngSpeedRanges(2,2)
                            headAngSpeedSample_loneWorm2{headAngSpeedSample_loneWorm2Ctr,1} = wormpath_xcoords_loneWorm;
                            headAngSpeedSample_loneWorm2{headAngSpeedSample_loneWorm2Ctr,2} = wormpath_ycoords_loneWorm;
                            sHeadAngleSpeedSample_loneWorm2{headAngSpeedSample_loneWorm2Ctr} = smoothHeadAngle_loneWorm;
                            if headAngSpeedSample_loneWorm2Ctr < size(headAngSpeedSample_loneWorm2,1)
                                headAngSpeedSample_loneWorm2Ctr = headAngSpeedSample_loneWorm2Ctr+1;
                            end
                            set(0,'CurrentFigure',HeadAngSpeedSample_loneWorm2Fig)
                            plot(wormpath_xcoords_loneWorm,wormpath_ycoords_loneWorm)
                        elseif sHeadAngSpeedRanges(3,1)<sHeadAngleChangeRate_loneWorm_thisFile(wormpathCtr,trajRunCtr) &...
                                sHeadAngleChangeRate_loneWorm_thisFile(wormpathCtr,trajRunCtr)<sHeadAngSpeedRanges(3,2)
                            headAngSpeedSample_loneWorm3{headAngSpeedSample_loneWorm3Ctr,1} = wormpath_xcoords_loneWorm;
                            headAngSpeedSample_loneWorm3{headAngSpeedSample_loneWorm3Ctr,2} = wormpath_ycoords_loneWorm;
                            sHeadAngleSpeedSample_loneWorm3{headAngSpeedSample_loneWorm3Ctr} = smoothHeadAngle_loneWorm;
                            if headAngSpeedSample_loneWorm3Ctr < size(headAngSpeedSample_loneWorm3,1)
                                headAngSpeedSample_loneWorm3Ctr = headAngSpeedSample_loneWorm3Ctr+1;
                            end
                            set(0,'CurrentFigure',HeadAngSpeedSample_loneWorm3Fig)
                            plot(wormpath_xcoords_loneWorm,wormpath_ycoords_loneWorm)
                        elseif sHeadAngSpeedRanges(4,1)<sHeadAngleChangeRate_loneWorm_thisFile(wormpathCtr,trajRunCtr) &...
                                sHeadAngleChangeRate_loneWorm_thisFile(wormpathCtr,trajRunCtr)<sHeadAngSpeedRanges(4,2)
                            headAngSpeedSample_loneWorm4{headAngSpeedSample_loneWorm4Ctr,1} = wormpath_xcoords_loneWorm;
                            headAngSpeedSample_loneWorm4{headAngSpeedSample_loneWorm4Ctr,2} = wormpath_ycoords_loneWorm;
                            sHeadAngleSpeedSample_loneWorm4{headAngSpeedSample_loneWorm4Ctr} = smoothHeadAngle_loneWorm;
                            if headAngSpeedSample_loneWorm4Ctr < size(headAngSpeedSample_loneWorm4,1)
                                headAngSpeedSample_loneWorm4Ctr = headAngSpeedSample_loneWorm4Ctr+1;
                            end
                            set(0,'CurrentFigure',HeadAngSpeedSample_loneWorm4Fig)
                            plot(wormpath_xcoords_loneWorm,wormpath_ycoords_loneWorm)
                        elseif sHeadAngSpeedRanges(5,1)<sHeadAngleChangeRate_loneWorm_thisFile(wormpathCtr,trajRunCtr) &...
                                sHeadAngleChangeRate_loneWorm_thisFile(wormpathCtr,trajRunCtr)<sHeadAngSpeedRanges(5,2)
                            headAngSpeedSample_loneWorm5{headAngSpeedSample_loneWorm5Ctr,1} = wormpath_xcoords_loneWorm;
                            headAngSpeedSample_loneWorm5{headAngSpeedSample_loneWorm5Ctr,2} = wormpath_ycoords_loneWorm;
                            sHeadAngleSpeedSample_loneWorm5{headAngSpeedSample_loneWorm5Ctr} = smoothHeadAngle_loneWorm;
                            if headAngSpeedSample_loneWorm5Ctr < size(headAngSpeedSample_loneWorm5,1)
                                headAngSpeedSample_loneWorm5Ctr = headAngSpeedSample_loneWorm5Ctr+1;
                            end
                            set(0,'CurrentFigure',HeadAngSpeedSample_loneWorm5Fig)
                            plot(wormpath_xcoords_loneWorm,wormpath_ycoords_loneWorm)
                        end
                    end
                end
            end
            headAngleChangeRate_leaveCluster_thisFile(isnan(headAngleChangeRate_leaveCluster_thisFile)) = [];
            headAngleChangeRate_loneWorm_thisFile(isnan(headAngleChangeRate_loneWorm_thisFile)) = [];
            sHeadAngleChangeRate_leaveCluster_thisFile(isnan(sHeadAngleChangeRate_leaveCluster_thisFile)) = [];
            sHeadAngleChangeRate_loneWorm_thisFile(isnan(sHeadAngleChangeRate_loneWorm_thisFile)) = [];
            % pool from different movies
            headAngleChangeRate_leaveCluster{fileCtr} = headAngleChangeRate_leaveCluster_thisFile;
            headAngleChangeRate_loneWorm{fileCtr} = headAngleChangeRate_loneWorm_thisFile;
            sHeadAngleChangeRate_leaveCluster{fileCtr} = sHeadAngleChangeRate_leaveCluster_thisFile;
            sHeadAngleChangeRate_loneWorm{fileCtr} = sHeadAngleChangeRate_loneWorm_thisFile;
            frameRunLengths_leaveCluster{fileCtr} = frameRunLengths_leaveCluster_thisFile;
            frameRunLengths_loneWorm{fileCtr} = frameRunLengths_loneWorm_thisFile;
        end
        % pool data from all files belonging to the same strain and worm density
        sHeadAngleChangeRate_leaveCluster = horzcat(sHeadAngleChangeRate_leaveCluster{:});
        sHeadAngleChangeRate_loneWorm = horzcat(sHeadAngleChangeRate_loneWorm{:});
        headAngleChangeRate_leaveCluster = horzcat(headAngleChangeRate_leaveCluster{:});
        headAngleChangeRate_loneWorm = horzcat(headAngleChangeRate_loneWorm{:});
        frameRunLengths_leaveCluster = horzcat(frameRunLengths_leaveCluster{:});
        frameRunLengths_loneWorm = horzcat(frameRunLengths_loneWorm{:});
        
        %% plot data, format, and export
        %
        set(0,'CurrentFigure',sHeadAngleChangeRateFig)
        histogram(sHeadAngleChangeRate_leaveCluster,'Normalization','pdf','DisplayStyle','stairs')
        hold on
        histogram(sHeadAngleChangeRate_loneWorm,'Normalization','pdf','DisplayStyle','stairs')
        leaveClusterLegend = strcat('leave cluster, n=',num2str(size(sHeadAngleChangeRate_leaveCluster,2)));
        loneWormLegend = strcat('lone worm, n=',num2str(size(sHeadAngleChangeRate_loneWorm,2)));
        legend(leaveClusterLegend, loneWormLegend)
        title([strains{strainCtr} '\_' wormnums{numCtr}],'FontWeight','normal')
        xlabel('head angle change rate (°/s)')
        ylabel('probability')
        %xlim([0 30])
        %ylim([0 0.09])
        set(sHeadAngleChangeRateFig,'PaperUnits','centimeters')
        figurename = ['figures/turns/headAngleChangeRate_' strains{strainCtr} '_' wormnums{numCtr} '_' phase '_data' num2str(dataset) '_' marker '_CL_2tDiv'];
        %savefig(sHeadAngleChangeRateFig,[figurename '.fig'])
        %exportfig(sHeadAngleChangeRateFig,[figurename '.eps'],exportOptions)
        %system(['epstopdf ' figurename '.eps']);
        %system(['rm ' figurename '.eps']);
        
        %
        set(0,'CurrentFigure',HeadAngSpeedSample_leaveCluster1Fig)
        title('HeadAngSpeedSample\_leaveCluster1Fig')
        set(0,'CurrentFigure',HeadAngSpeedSample_leaveCluster2Fig)
        title('HeadAngSpeedSample\_leaveCluster2Fig')
        set(0,'CurrentFigure',HeadAngSpeedSample_leaveCluster3Fig)
        title('HeadAngSpeedSample\_leaveCluster3Fig')
        set(0,'CurrentFigure',HeadAngSpeedSample_leaveCluster4Fig)
        title('HeadAngSpeedSample\_leaveCluster4Fig')
        set(0,'CurrentFigure',HeadAngSpeedSample_leaveCluster5Fig)
        title('HeadAngSpeedSample\_leaveCluster5Fig')
        set(0,'CurrentFigure',HeadAngSpeedSample_loneWorm1Fig)
        title('HeadAngSpeedSample\_loneWorm1Fig')
        set(0,'CurrentFigure',HeadAngSpeedSample_loneWorm1Fig)
        title('HeadAngSpeedSample\_loneWorm1Fig')
        set(0,'CurrentFigure',HeadAngSpeedSample_loneWorm2Fig)
        title('HeadAngSpeedSample\_loneWorm2Fig')
        set(0,'CurrentFigure',HeadAngSpeedSample_loneWorm3Fig)
        title('HeadAngSpeedSample\_loneWorm3Fig')
        set(0,'CurrentFigure',HeadAngSpeedSample_loneWorm4Fig)
        title('HeadAngSpeedSample\_loneWorm4Fig')
        set(0,'CurrentFigure',HeadAngSpeedSample_loneWorm5Fig)
        title('HeadAngSpeedSample\_loneWorm5Fig')
        %save(['figures/turns/results/sHeadAngleSpeedSample_leaveCluster1_'...
        %    strains{strainCtr} '_' wormnums{numCtr} '_' phase '_data' num2str(dataset) '_' marker '.mat'],...
        %    'sHeadAngleSpeedSample_leaveCluster1')
    end
end