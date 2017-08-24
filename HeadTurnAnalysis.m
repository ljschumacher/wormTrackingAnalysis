% Script plots three head angular features: 
% 1. total head angle change over the full trajectory, 
% 2. normalised head angle change over full trajectory (normalised
% by path length), and 
% 3. head angular speed. 
% Script compares leaveCluster and loneWorm, and can be modified to include further worm categories. 
% Script also contains an option to visualise sample trajectories that give rise
% to values that fall within pre-specified value ranges. 

clear
close all

%% set parameters
phase = 'fullMovie'; % 'fullMovie', 'joining', or 'sweeping'.
dataset = 2; % 1 or 2
marker = 'bodywall'; % 'pharynx' or 'bodywall'
strains = {'npr1','N2'}; % {'npr1','N2'}
wormnums = {'40'};% {'40'};
wormcats = {'leaveCluster','loneWorm'}; %'leaveCluster','loneWorm'
smoothing = false;
postExitDuration = 5; % duration (in seconds) after a worm exits a cluster to be included in the leave cluster analysis
minTrajDuration = 1; % duration (in seconds) of minimum traj length
maxTrajDuration = 5;  % duration (in seconds) of maximum traj length % may set to 1.5 to truncate loneWorm traj to match those of leaveCluster traj length
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
saveResults = false;
visualiseSampleTraj = true; % true or false

if visualiseSampleTraj == true
    numSampleTraj = 5; % number of sample trajectories to be plotted for each angular speed range
    featureToSample = 'headAngTotal'; % 'headAngTotal','headAngNorm', or 'headAngSpeed'
    if strcmp(featureToSample,'headAngTotal') | strcmp(featureToSample,'headAngSpeed')
        headAngRanges = [0, 0.25; pi/2-0.25, pi/2+0.25; pi-0.25, pi+0.25; 3/2*pi-0.25, 3/2*pi+0.25; 2*pi-0.25, 2*pi];
    elseif strcmp(featureToSample,'headAngNorm')
        headAngRanges = [0 0.01; 0.01 0.02; 0.02 0.03; 0.03 0.05; 0.05 1];
    else
        warning('Wrong feature selected for trajectory visualisation')
    end
end

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
            headAngTotal.(wormcats{wormcatCtr}) = cell(numFiles,1);
            headAngNorm.(wormcats{wormcatCtr}) = cell(numFiles,1);
            headAngSpeed.(wormcats{wormcatCtr}) = cell(numFiles,1);
            frameRunLengths.(wormcats{wormcatCtr}) = cell(numFiles,1);
        end
        
        if visualiseSampleTraj
            % save up to 500 sets of xy coordinates for paths that fall within a certain angular speed range to plot sample trajectories
            for wormcatCtr = 1:length(wormcats)
                headAngSampleTraj.(wormcats{wormcatCtr}) = cell(500,2,size(headAngRanges,1));
                headAngSampleCtr.(wormcats{wormcatCtr}) = ones(size(headAngRanges,1),1);
            end
        end
        
        %% go through individual movies
        for fileCtr = 1:numFiles
            fileCtr
            %% load data
            filename = filenames{fileCtr};
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
            if strcmp(marker,'bodyWall')
                worm_xcoords = worm_xcoords(:,1:8); % restrict to head nodes only (8 out of 49 nodes)
                worm_ycoords = worm_ycoords(:,1:8);
            end
            uniqueWorms = unique(trajData.worm_index_joined);
            
            % initialise
            for wormcatCtr = 1:length(wormcats)
                % assume each worm has up to 100 leave cluster traj. Checks for this later and warns if not enough
                headAngTotal.(wormcats{wormcatCtr}){fileCtr} = NaN(numel(uniqueWorms),100);
                headAngNorm.(wormcats{wormcatCtr}){fileCtr} = NaN(numel(uniqueWorms),100);
                headAngSpeed.(wormcats{wormcatCtr}){fileCtr} = NaN(numel(uniqueWorms),100);
                % create variable to keep track of trajectory lengths
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
                        continuousFramesMinLengthLogInd = cellfun(@numel,continuousFrameRuns)>=round(minTrajDuration*frameRate);
                        continuousFrameRuns = continuousFrameRuns(continuousFramesMinLengthLogInd);
                        
                        % go through each traj that fits the min length criteria to obtain xy coordinates
                        if size(continuousFrameRuns,2)>100
                            warning('more trajectories present than allocated space to hold values for')
                        end
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
                            maxTrajFrameNum = round(maxTrajDuration*frameRate);
                            if size(wormtraj_xcoords,1)> maxTrajFrameNum
                                if strcmp(wormcats{wormcatCtr},'leaveCluster')
                                    % always start leave cluster traj from the moment the worm exits cluster
                                    firstFrame = 1;
                                elseif strcmp(wormcats{wormcatCtr},'loneWorm')
                                    % randomly sample the start of lone worm traj
                                    firstFrame = randi(size(wormtraj_xcoords,1)-maxTrajFrameNum,1);
                                end
                                % truncate the traj at maximum length
                                lastFrame = firstFrame + maxTrajFrameNum;
                                wormtraj_xcoords = wormtraj_xcoords(firstFrame:lastFrame,:);
                                wormtraj_ycoords = wormtraj_ycoords(firstFrame:lastFrame,:);
                            end
                            
                            % calculate head angle changes per frame
                            [headAngleDiff, framesElapsed] = getHeadAngleDiff(wormtraj_xcoords,wormtraj_ycoords, marker, smoothing, frameRate);
                            % calculate total head turn over full trajectory
                            headAngTotal.(wormcats{wormcatCtr}){fileCtr}(wormCtr,trajCtr) =...
                                abs(nansum(headAngleDiff));
                            % normalise head turn over path length
                            wormtraj_xcoords = mean(wormtraj_xcoords,2)*pixelsize; % turn pixels into microns
                            wormtraj_ycoords = mean(wormtraj_ycoords,2)*pixelsize;
                            xdiff = wormtraj_xcoords(2:end) - wormtraj_xcoords(1:end-1);
                            ydiff = wormtraj_ycoords(2:end) - wormtraj_ycoords(1:end-1);
                            pathLength = sum(sqrt(xdiff.^2+ydiff.^2)); % calculate path length
                            headAngNorm.(wormcats{wormcatCtr}){fileCtr}(wormCtr,trajCtr) =...
                                headAngTotal.(wormcats{wormcatCtr}){fileCtr}(wormCtr,trajCtr)/pathLength;
                            % calculate head angular speed
                            headAngSpeed.(wormcats{wormcatCtr}){fileCtr}(wormCtr,trajCtr) = ...
                                headAngTotal.(wormcats{wormcatCtr}){fileCtr}(wormCtr,trajCtr)/framesElapsed*frameRate;
                            
                            % optional: save xy coordinates for paths that fall within certain head angle measurement ranges
                            if visualiseSampleTraj
                                % loop through each range to see which one it falls within
                                for rangeCtr = 1:size(headAngRanges,1)
                                    % if a headAngSpeed value falls in between the limits of specified ranges
                                    if headAngRanges(rangeCtr,1)<headAngTotal.(wormcats{wormcatCtr}){fileCtr}(wormCtr,trajCtr) & ...
                                            headAngTotal.(wormcats{wormcatCtr}){fileCtr}(wormCtr,trajCtr)<headAngRanges(rangeCtr,2)
                                        % save xy coordinates (in microns)
                                        headAngSampleTraj.(wormcats{wormcatCtr}){headAngSampleCtr.(wormcats{wormcatCtr})(rangeCtr),1,rangeCtr} = wormtraj_xcoords;
                                        headAngSampleTraj.(wormcats{wormcatCtr}){headAngSampleCtr.(wormcats{wormcatCtr})(rangeCtr),2,rangeCtr} = wormtraj_ycoords;
                                        % update the traj counter for that range (until up to 500 preallocated trajectory spaces)
                                        if headAngSampleCtr.(wormcats{wormcatCtr})(rangeCtr) < size(headAngSampleTraj.(wormcats{wormcatCtr}),1)
                                            headAngSampleCtr.(wormcats{wormcatCtr})(rangeCtr) = headAngSampleCtr.(wormcats{wormcatCtr})(rangeCtr)+1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            for wormcatCtr = 1:length(wormcats)
                % remove empty entries
                headAngTotal.(wormcats{wormcatCtr}){fileCtr}(isnan(headAngTotal.(wormcats{wormcatCtr}){fileCtr}))=[];
                headAngNorm.(wormcats{wormcatCtr}){fileCtr}(isnan(headAngNorm.(wormcats{wormcatCtr}){fileCtr}))=[];
                headAngSpeed.(wormcats{wormcatCtr}){fileCtr}(isnan(headAngSpeed.(wormcats{wormcatCtr}){fileCtr}))=[];
            end
        end
        % pool data from all files belonging to the same strain and worm density
        for wormcatCtr = 1:length(wormcats)
            headAngTotalPool.(wormcats{wormcatCtr}) = horzcat(headAngTotal.(wormcats{wormcatCtr}){:});
            headAngNormPool.(wormcats{wormcatCtr}) = horzcat(headAngNorm.(wormcats{wormcatCtr}){:});
            headAngSpeedPool.(wormcats{wormcatCtr}) = horzcat(headAngSpeed.(wormcats{wormcatCtr}){:});
            frameRunLengths.(wormcats{wormcatCtr}) = horzcat(frameRunLengths.(wormcats{wormcatCtr}){:});
        end
        
        %% plot data, format, and export
        % save head angle values
        if saveResults
            save(['figures/turns/results/headAngTotal_' strain '_' wormnum '_' phase '_data' num2str(dataset) '_' marker '0.mat'],'headAngTotal')
            save(['figures/turns/results/headAngNorm_' strain '_' wormnum '_' phase '_data' num2str(dataset) '_' marker '0.mat'],'headAngNorm')
            save(['figures/turns/results/headAngSpeed_' strain '_' wormnum '_' phase '_data' num2str(dataset) '_' marker '0.mat'],'headAngSpeed')
            save(['figures/turns/results/frameRunLengths_' strain '_' wormnum '_' phase '_data' num2str(dataset) '_' marker '0.mat'],'frameRunLengths')
        end
        
        % plot total head angle change
        headAngTotalFig = figure; hold on
        for wormcatCtr = 1:length(wormcats)
            histogram(headAngTotalPool.(wormcats{wormcatCtr}),'Normalization','pdf','DisplayStyle','stairs')
            Legend{wormcatCtr} = strcat(wormcats{wormcatCtr}, ',n=', num2str(size(headAngTotalPool.(wormcats{wormcatCtr}),2)));
        end
        legend(Legend)
        title([strains{strainCtr} '\_' wormnum],'FontWeight','normal')
        xlabel('Total head angle change (radian)')
        ylabel('Probability')
        xlim([0 7])
        %ylim([0 2])
        set(headAngTotalFig,'PaperUnits','centimeters')
        figurename = ['figures/turns/headAngTotal_' strain '_' wormnum '_' phase '_data' num2str(dataset) '_' marker];
        if saveResults
            load('exportOptions.mat')
            exportfig(headAngTotalFig,[figurename '0.eps'],exportOptions)
            system(['epstopdf ' figurename '.eps']);
            system(['rm ' figurename '.eps']);
        end
        
        % plot normalised head angle change
        headAngNormFig = figure; hold on
        for wormcatCtr = 1:length(wormcats)
            histogram(headAngNormPool.(wormcats{wormcatCtr}),'Normalization','pdf','DisplayStyle','stairs')
            Legend{wormcatCtr} = strcat(wormcats{wormcatCtr}, ',n=', num2str(size(headAngNormPool.(wormcats{wormcatCtr}),2)));
        end
        legend(Legend)
        title([strains{strainCtr} '\_' wormnum],'FontWeight','normal')
        xlabel('Normalised head angle change (radian/micron)')
        ylabel('Probability')
        if strcmp(marker,'pharynx')
            xlim([0 0.1])
        elseif strcmp(marker,'bodywall')
            xlim([0 0.3])
        end
        %ylim([0 2])
        set(headAngNormFig,'PaperUnits','centimeters')
        figurename = ['figures/turns/headAngNorm_' strain '_' wormnum '_' phase '_data' num2str(dataset) '_' marker];
        if saveResults
            load('exportOptions.mat')
            exportfig(headAngNormFig,[figurename '0.eps'],exportOptions)
            system(['epstopdf ' figurename '.eps']);
            system(['rm ' figurename '.eps']);
        end
        
        % plot head angle speed
        headAngSpeedFig = figure; hold on
        for wormcatCtr = 1:length(wormcats)
            histogram(headAngSpeedPool.(wormcats{wormcatCtr}),'Normalization','pdf','DisplayStyle','stairs')
            Legend{wormcatCtr} = strcat(wormcats{wormcatCtr}, ',n=', num2str(size(headAngSpeedPool.(wormcats{wormcatCtr}),2)));
        end
        legend(Legend)
        title([strain '\_' wormnum],'FontWeight','normal')
        xlabel('Head angular speed (radian/s)')
        ylabel('Probability')
        xlim([0 7])
        %ylim([0 2])
        set(headAngSpeedFig,'PaperUnits','centimeters')
        figurename = ['figures/turns/headAngSpeed_' strain '_' wormnum '_' phase '_data' num2str(dataset) '_' marker];
        if saveResults
            load('exportOptions.mat')
            exportfig(headAngSpeedFig,[figurename '0.eps'],exportOptions)
            system(['epstopdf ' figurename '.eps']);
            system(['rm ' figurename '.eps']);
        end
        
        if visualiseSampleTraj
            % save trajectories
            if saveResults
            save(['figures/turns/results/headAngSampleTraj_' strain '_' wormnum '_' phase '_data' num2str(dataset) '_' marker '.mat'],'headAngSampleTraj')
            end
            % plot trajectories
            plotSampleHeadAngTraj(headAngSampleTraj,headAngRanges,featureToSample,numSampleTraj,wormcats,strain,wormnum,marker,phase, dataset, saveResults)
        end
    end
end