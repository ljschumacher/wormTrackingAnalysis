clear
close all

% Script plots midbody speed for worms entering or leaving a cluster using the manually joined red trajectories, with
% the point of entry/exit aligned.

%% set parameters
phase = 'joining'; % 'fullMovie', 'joining', or 'sweeping'.
strains = {'npr1'}; % {'npr1','N2'}
wormnums = {'40'};% {'40'};
preExitDuration = 20; % duration (in seconds) before a worm exits a cluster to be included in the leave cluster analysis
postExitDuration = 20; % duration (in seconds) after a worm exits a cluster to be included in the leave cluster analysis
smoothWindow = 9; % number of frames to smooth over
saveResults = true;

useManualEvents = true; % manual events only available for joining phase
if useManualEvents
    manualEventMaxDuration = 400; % max number of frames that contains the beginning and end of an annotated event
    maxTrajPerGraph = 10;
    
else
    minInOutClusterFrameNum = 5; % unless using manually labelled events, then enter the min number of frames for in/out cluster status to be sustained for an event to be considered
    enforceInClusterAfterEntryBeforeExit = true;
end

pixelsize = 100/19.5; % 100 microns are 19.5 pixels
maxBlobSize_r = 2.5e5;
minSkelLength_r = 850;
maxSkelLength_r = 1500;
minNeighbrDist = 2000;
inClusterNeighbourNum = 3;

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',30,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',20,...
    'LineWidth',3);

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        
        %% load file list and annotations
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_r_list.xlsx'],1,'A1:E15','basic');
        numFiles = length(filenames);
        if useManualEvents
            if strcmp(strains{strainCtr},'npr1')
                [~,~,annotations] = xlsread(['/data2/shared/data/twoColour/MaskedVideos/entryExitEvents_' strains{strainCtr} '_' wormnum '_' phase '.xlsx'],1,'A2:I80','basic');
            end
            % count the total number of entry and exit events
            totalEntry = 0;
            totalExit = 0;
            for eventCtr = 1:size(annotations,1)
                if strcmp(annotations{eventCtr,5},'enter')
                    totalEntry = totalEntry + 1;
                    entry(totalEntry) = eventCtr;
                elseif strcmp(annotations{eventCtr,5},'exit')
                    totalExit = totalExit+1;
                end
            end
            % pre-allocate matrix to write speeds for each event
            frameRate = 9; % should be the case for most if not all files..
            entrySpeeds = NaN(totalEntry,frameRate*(preExitDuration+postExitDuration)+1+manualEventMaxDuration); % assuming no more than 400 frames for the actual entry/exit event
            exitSpeeds = NaN(totalExit,frameRate*(preExitDuration+postExitDuration)+1+manualEventMaxDuration);
        end
        
        % initialise counters
        entryCtr = 1;
        exitCtr = 1;
        %% go through individual movies
        for fileCtr = 1:numFiles
            %% load data
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton'); % in pixels
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            features = h5read(strrep(filename,'skeletons','feat_manual'),'/features_timeseries');
            
            %% generate time series (x), setting time to be zero at the point of entry start
            if useManualEvents
                timeSeries.entry = [-preExitDuration*frameRate:postExitDuration*frameRate+manualEventMaxDuration];
                timeSeries.exit = [-preExitDuration*frameRate-manualEventMaxDuration:postExitDuration*frameRate];
            else
                timeSeries = [-preExitDuration*frameRate:postExitDuration*frameRate];
            end
            
            %% calculate midbody signed speed (from reversalAnalysisBodyWall.m)
            midbodyIndcs = 19:33;
            % centroids of midbody skeleton
            midbody_x = mean(squeeze(skelData(1,midbodyIndcs,:)))*pixelsize;
            midbody_y = mean(squeeze(skelData(2,midbodyIndcs,:)))*pixelsize;
            % change in centroid position over time (issue of worm shifts?)
            dmidbody_xdt = gradient(midbody_x)*frameRate;
            dmidbody_ydt = gradient(midbody_y)*frameRate;
            % midbody speed and velocity
            dFramedt = gradient(double(trajData.frame_number))';
            midbodySpeed = sqrt(dmidbody_xdt.^2 + dmidbody_ydt.^2)./dFramedt;
            midbodyVelocity = [dmidbody_xdt; dmidbody_ydt]./dFramedt;
            % direction of segments pointing along midbody
            [~, dmidbody_yds] = gradient(squeeze(skelData(2,midbodyIndcs,:)),-1);
            [~, dmidbody_xds] = gradient(squeeze(skelData(1,midbodyIndcs,:)),-1);
            % sign speed based on relative orientation of velocity to midbody
            midbodySpeedSigned = getSignedSpeed(midbodyVelocity,[mean(dmidbody_xds); mean(dmidbody_yds)]);
            % ignore first and last frames of each worm's track
            wormChangeIndcs = gradient(double(trajData.worm_index_manual))~=0;
            midbodySpeedSigned(wormChangeIndcs)=NaN;
            
            %% filter worms by various criteria
            % filter red by manually joined traj
            trajData.filtered = ismember(trajData.worm_index_manual,int32(features.worm_index));
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
            % apply phase restriction
            [firstPhaseFrame, lastPhaseFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr);
            phaseFrameLogInd = trajData.frame_number <= lastPhaseFrame & trajData.frame_number >= firstPhaseFrame;
            trajData.filtered(~phaseFrameLogInd) = false;
            if ~useManualEvents
                % get logical indices for worms that have entered or left a cluster
                min_neighbr_dist = h5read(filename,'/min_neighbr_dist');
                num_close_neighbrs = h5read(filename,'/num_close_neighbrs');
                trajData = h5read(filename,'/trajectories_data');
                frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
                inClusterLogInd = num_close_neighbrs>=inClusterNeighbourNum;
                enterClusterStartLogInd = vertcat(false,~inClusterLogInd(1:end-1)&inClusterLogInd(2:end));
                exitClusterStartLogInd = vertcat(false,inClusterLogInd(1:end-1)&~inClusterLogInd(2:end));
            end
            
            %% obtain speed over time course for cluster entry events
            
            %% case with manually labelled events
            if useManualEvents
                % go through each annotated event
                for eventCtr = 1:size(annotations,1)
                    % find entry events with matching recording file name to the file currently under consideration
                    recordingNumber = annotations{eventCtr,1};
                    recordingNumber = recordingNumber(2:end); % remove the 'r' before the number
                    if contains(filename, recordingNumber) & strcmp(annotations{eventCtr,5},'enter')
                        % collect information on worm index and entry frames
                        wormIndex = annotations{eventCtr,2};
                        thisEntryStartFrame = annotations{eventCtr,7}; % get the annotated entry start frame
                        thisEntryEndFrame= annotations{eventCtr,8}; % get the annotated entry finish frame
                        if thisEntryEndFrame > lastPhaseFrame % trim hand annotation in case event end goes beyond phase end
                            thisEntryEndFrame = lastPhaseFrame;
                        end
                        if thisEntryEndFrame-thisEntryStartFrame +1 > manualEventMaxDuration
                            warning('event duration exceeds allocated maximum duration, increase value for manualEventMaxDuration variable')
                        end
                        % extend for a specified duration before and after the entry point
                        thisEntryXStartFrame = thisEntryStartFrame - preExitDuration*frameRate;
                        if thisEntryXStartFrame <= firstPhaseFrame
                            beforeStartFrameNum = firstPhaseFrame-thisEntryXStartFrame; % take note of omitted frames for alignment purposes
                            thisEntryXStartFrame = firstPhaseFrame; % exclude frames before the start of the specified phase
                        end
                        thisEntryXEndFrame = thisEntryEndFrame + postExitDuration*frameRate;
                        if  thisEntryXEndFrame > lastPhaseFrame
                            afterEndFrameNum =  thisEntryXEndFrame -lastPhaseFrame; % take note of omitted frames for alignment purposes
                            thisEntryXEndFrame = lastPhaseFrame; % exclude frames beyond the end of the specified phase
                        end
                        % get aligned list of frames for the event
                        thisEntrySpeeds = NaN(1,frameRate*(preExitDuration+postExitDuration)+1+manualEventMaxDuration);
                        thisEntryXFrames = NaN(1,frameRate*(preExitDuration+postExitDuration)+1+manualEventMaxDuration);
                        if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                            thisEntryXFrames(beforeStartFrameNum+1:end) = thisEntryXStartFrame:thisEntryXEndFrame;
                        elseif exist('afterEndFrameNum','var')
                            thisEntryXFrames(1:(end-afterEndFrameNum)) = thisEntryXStartFrame:thisEntryXEndFrame;
                        else
                            thisEntryXFrames = thisEntryXStartFrame:thisEntryXEndFrame;
                        end
                        % go through each frame
                        for frameCtr = 1:length(thisEntryXFrames)
                            frameNumber = thisEntryXFrames(frameCtr);
                            if ~isnan(frameNumber)
                                wormFrameLogInd = trajData.filtered & ...
                                trajData.worm_index_manual == wormIndex & trajData.frame_number == frameNumber;
                                if nnz(wormFrameLogInd)~=0
                                    assert(nnz(wormFrameLogInd) ==1);
                                    thisEntrySpeeds(frameCtr) = midbodySpeedSigned(wormFrameLogInd);
                                end
                            end
                        end
                        % write speeds
                        entrySpeeds(entryCtr,1:length(thisEntrySpeeds)) = thisEntrySpeeds;
                        % clear variables
                        clear beforeStartFrameNum
                        clear afterEndFrameNum
                        % write event number to legend
                        entryLegend{entryCtr} = num2str(eventCtr);
                        % update counter
                        entryCtr = entryCtr +1;
                        assert(entryCtr<=totalEntry+1);
                    end
                end
                
                %% without manual event annotation, we need to identify our own events
            else
                % get indices for cluster entry points
                enterClusterPoint = find(enterClusterStartLogInd);
                entrySpeeds = NaN(numel(enterClusterPoint),preExitDuration*frameRate +1+ postExitDuration*frameRate);
                % loop through each entry
                for entryCtr = 1:numel(enterClusterPoint)
                    thisEntryIdx = enterClusterPoint(entryCtr);
                    % get the worm index
                    wormIndex = trajData.worm_index_manual(thisEntryIdx);
                    % check that the same worm stays in cluster for minumum number of frames
                    if thisEntryIdx+minInOutClusterFrameNum <= length(trajData.frame_number) &&...
                            nnz(inClusterLogInd(thisEntryIdx:(thisEntryIdx+minInOutClusterFrameNum))) == 1+minInOutClusterFrameNum...
                            && nnz(trajData.worm_index_manual(thisEntryIdx:(thisEntryIdx+minInOutClusterFrameNum)) == wormIndex) == 1+minInOutClusterFrameNum
                        % expand for a specified number of entries before and after the entry point
                        thisEntryStartIdx = thisEntryIdx-preExitDuration*frameRate;
                        if thisEntryStartIdx <= 0
                            beforeStartFrameNum = 1-thisEntryStartIdx; % take note of omitted frames for alignment purposes
                            thisEntryStartIdx = 1; % exclude entries below index 0
                        end
                        thisEntryEndIdx = thisEntryIdx + postExitDuration*frameRate;
                        if thisEntryEndIdx > length(trajData.frame_number)
                            afterEndFrameNum = thisEntryEndIdx - length(trajData.frame_number); % take note of omitted frames for alignment purposes
                            thisEntryEndIdx = length(trajData.frame_number); % exclude entries above highest index
                        end
                        % generate logical index for this entry
                        thisEntryLogInd = false(1,length(midbodySpeedSigned));
                        thisEntryLogInd(thisEntryStartIdx:thisEntryEndIdx) = true;
                        % get midbody signed speeds (y) for this entry
                        thisEntryXSpeeds = NaN(1,size(entrySpeeds,2));
                        if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                            thisEntryXSpeeds(beforeStartFrameNum+1:end) = midbodySpeedSigned(thisEntryLogInd);
                        elseif exist('afterEndFrameNum','var')
                            thisEntryXSpeeds(1:(end-afterEndFrameNum)) = midbodySpeedSigned(thisEntryLogInd);
                        else
                            thisEntryXSpeeds = midbodySpeedSigned(thisEntryLogInd);
                        end
                        % exclude entries representing a different worm
                        thisEntryXSpeeds(trajData.worm_index_manual(thisEntryLogInd) ~= wormIndex) = NaN;
                        % optional: excludes entries representing non-inCluster worms post-entry (inclusive of entry point)
                        if enforceInClusterAfterEntryBeforeExit
                            % create speed vector the same length as the full index
                            thisEntrySpeedsFullLength = NaN(1,length(midbodySpeedSigned));
                            if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                                thisEntrySpeedsFullLength(thisEntryStartIdx:thisEntryEndIdx) = thisEntryXSpeeds(beforeStartFrameNum+1:end);
                            elseif exist('afterEndFrameNum','var')
                                thisEntrySpeedsFullLength(thisEntryStartIdx:thisEntryEndIdx) = thisEntryXSpeeds(1:(end-afterEndFrameNum));
                            else
                                thisEntrySpeedsFullLength(thisEntryStartIdx:thisEntryEndIdx) = thisEntryXSpeeds;
                            end
                            % generate logical index for post-entry indices
                            thisEntryPostEntryLogInd = false(1,length(midbodySpeedSigned));
                            thisEntryPostEntryLogInd(thisEntryIdx:thisEntryEndIdx) = true;
                            % eliminate speeds for non-inCluster worms post-entry
                            thisEntrySpeedsFullLength(thisEntryPostEntryLogInd & ~ inClusterLogInd') = NaN;
                            % shorten full length speeds back to the window of interest
                            if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                                thisEntryXSpeeds(beforeStartFrameNum+1:end) = thisEntrySpeedsFullLength(thisEntryStartIdx:thisEntryEndIdx);
                            elseif exist('afterEndFrameNum','var')
                                thisEntryXSpeeds(1:(end-afterEndFrameNum)) = thisEntrySpeedsFullLength(thisEntryStartIdx:thisEntryEndIdx);
                            else
                                thisEntryXSpeeds = thisEntrySpeedsFullLength(thisEntryStartIdx:thisEntryEndIdx);
                            end
                        end
                        % add data to speed matrix
                        entrySpeeds(entryCtr,:) = thisEntryXSpeeds;
                    end
                    % clear variables
                    clear beforeStartFrameNum
                    clear afterEndFrameNum
                end
            end
            
            % set maximum speed and remove 0 speed
            entrySpeeds(abs(entrySpeeds)>1500) = NaN;
            entrySpeeds(entrySpeeds==0) = NaN;
            % smooth speeds
            smoothEntrySpeeds = smoothdata(entrySpeeds,2,'movmean',smoothWindow,'includenan');
            % set maximum speed
            smoothEntrySpeeds(abs(smoothEntrySpeeds)>1500) = NaN;
            
            %% obtain speed over time course for cluster exit events
            
            %% case with manually labelled events
            if useManualEvents
                % go through each annotated event
                for eventCtr = 1:size(annotations,1)
                    % find exit events with matching recording file name to the file currently under consideration
                    recordingNumber = annotations{eventCtr,1};
                    recordingNumber = recordingNumber(2:end); % remove the 'r' before the number
                    if contains(filename, recordingNumber) & strcmp(annotations{eventCtr,5},'exit')
                        % collect information on worm index and exit frames
                        wormIndex = annotations{eventCtr,2};
                        thisExitStartFrame = annotations{eventCtr,7}; % get the annotated entry start frame
                        thisExitEndFrame= annotations{eventCtr,8}; % get the annotated entry finish frame
                        if thisExitEndFrame > lastPhaseFrame % trim hand annotation in case event end goes beyond phase end
                            thisExitEndFrame = lastPhaseFrame;
                        end
                        if thisExitEndFrame-thisExitStartFrame +1 > manualEventMaxDuration
                            warning('event duration exceeds allocated maximum duration, increase value for manualEventMaxDuration variable')
                        end
                        % extend for a specified duration before and after the exit point
                        thisExitXStartFrame = thisExitStartFrame - preExitDuration*frameRate;
                        if thisExitXStartFrame <= firstPhaseFrame
                            beforeStartFrameNum = firstPhaseFrame-thisExitXStartFrame; % take note of omitted frames for alignment purposes
                            thisExitXStartFrame = firstPhaseFrame; % exclude frames before the start of the specified phase
                        end
                        thisExitXEndFrame = thisExitEndFrame + postExitDuration*frameRate;
                        if  thisExitXEndFrame > lastPhaseFrame
                            afterEndFrameNum =  thisExitXEndFrame -lastPhaseFrame; % take note of omitted frames for alignment purposes
                            thisExitXEndFrame = lastPhaseFrame; % exclude frames beyond the end of the specified phase
                        end
                        % get aligned list of frames for the event
                        exitNumFrames = thisExitEndFrame - thisExitStartFrame;
                        startFiller = manualEventMaxDuration - exitNumFrames ; % number of empty frames to add to the start of speed vector to keep alignment for end of exit
                        thisExitSpeeds = NaN(1,frameRate*(preExitDuration+postExitDuration)+1+manualEventMaxDuration);
                        thisExitXFrames = NaN(1,frameRate*(preExitDuration+postExitDuration)+1+manualEventMaxDuration);
                        if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                            thisExitXFrames(beforeStartFrameNum+startFiller+1:end) = thisExitXStartFrame:thisExitXEndFrame;
                        elseif exist('afterEndFrameNum','var')
                            thisExitXFrames(startFiller+1:(end-afterEndFrameNum)) = thisExitXStartFrame:thisExitXEndFrame;
                        else
                            thisExitXFrames(startFiller+1:end) = thisExitXStartFrame:thisExitXEndFrame;
                        end
                        % go through each frame
                        for frameCtr = 1:length(thisExitXFrames)
                            frameNumber = thisExitXFrames(frameCtr);
                            if ~isnan(frameNumber)
                                wormFrameLogInd = trajData.filtered &...
                                trajData.worm_index_manual == wormIndex & trajData.frame_number == frameNumber;
                                if nnz(wormFrameLogInd)~=0
                                    assert(nnz(wormFrameLogInd) ==1);
                                    thisExitSpeeds(frameCtr) = midbodySpeedSigned(wormFrameLogInd);
                                end
                            end
                        end
                        % write speeds
                        exitSpeeds(exitCtr,1:length(thisExitSpeeds)) = thisExitSpeeds;
                        % clear variables
                        clear beforeStartFrameNum
                        clear afterEndFrameNum
                        % write event number to legend
                        exitLegend{exitCtr} = num2str(eventCtr);
                        % update exit counter
                        exitCtr = exitCtr + 1;
                        assert(exitCtr<=totalExit+1);
                    end
                end
                
                %% without manual event annotation, we need to identify our own events
            else
                % get indices for cluster exit points
                exitClusterPoint = find(exitClusterStartLogInd);
                exitSpeeds = NaN(numel(exitClusterPoint),preExitDuration*frameRate +1+ postExitDuration*frameRate);
                % loop through each exit
                for exitCtr = 1:numel(exitClusterPoint)
                    thisExitIdx = exitClusterPoint(exitCtr);
                    % get the worm index
                    wormIndex = trajData.worm_index_manual(thisExitIdx);
                    % check that the same worm stays out of cluster for minumum number of frames
                    if thisExitIdx+minInOutClusterFrameNum <= length(trajData.frame_number) &&...
                            nnz(inClusterLogInd(thisExitIdx:(thisExitIdx+minInOutClusterFrameNum))) == 0 ...
                            && nnz(trajData.worm_index_manual(thisExitIdx:(thisExitIdx+minInOutClusterFrameNum)) == wormIndex) == 1+minInOutClusterFrameNum
                        % expand for a specified number of entries before and after the exit point
                        thisExitStartIdx = thisExitIdx-preExitDuration*frameRate;
                        if thisExitStartIdx <= 0
                            beforeStartFrameNum = 1-thisExitStartIdx; % take note of omitted frames for alignment purposes
                            thisExitStartIdx = 1; % exclude entries below index 0
                        end
                        thisExitEndIdx = thisExitIdx + postExitDuration*frameRate;
                        if thisExitEndIdx > length(trajData.frame_number)
                            afterEndFrameNum = thisExitEndIdx - length(trajData.frame_number); % take note of omitted frames for alignment purposes
                            thisExitEndIdx = length(trajData.frame_number); % exclude entries above highest index
                        end
                        % generate logical index for this exit
                        thisExitLogInd = false(1,length(midbodySpeedSigned));
                        thisExitLogInd(thisExitStartIdx:thisExitEndIdx) = true;
                        % get midbody signed speeds (y) for this exit
                        thisExitSpeeds = NaN(1,size(exitSpeeds,2));
                        if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                            thisExitSpeeds(beforeStartFrameNum+1:end) = midbodySpeedSigned(thisExitLogInd);
                        elseif exist('afterEndFrameNum','var')
                            thisExitSpeeds(1:(end-afterEndFrameNum)) = midbodySpeedSigned(thisExitLogInd);
                        else
                            thisExitSpeeds = midbodySpeedSigned(thisExitLogInd);
                        end
                        % exclude entries representing a different worm
                        thisExitSpeeds(trajData.worm_index_manual(thisExitLogInd) ~= wormIndex) = NaN;
                        % optional: excludes entries representing non-inCluster worms post-entry (inclusive of entry point)
                        if enforceInClusterAfterEntryBeforeExit
                            % create speed vector the same length as the full index
                            thisExitSpeedsFullLength = NaN(1,length(midbodySpeedSigned));
                            if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                                thisExitSpeedsFullLength(thisExitStartIdx:thisExitEndIdx) = thisExitSpeeds(beforeStartFrameNum+1:end);
                            elseif exist('afterEndFrameNum','var')
                                thisExitSpeedsFullLength(thisExitStartIdx:thisExitEndIdx) = thisExitSpeeds(1:(end-afterEndFrameNum));
                            else
                                thisExitSpeedsFullLength(thisExitStartIdx:thisExitEndIdx) = thisExitSpeeds;
                            end
                            % generate logical index for post-entry indices
                            thisExitPreExitLogInd = false(1,length(midbodySpeedSigned));
                            thisExitPreExitLogInd(thisExitStartIdx:(thisExitIdx-1)) = true;
                            % eliminate speeds for non-inCluster worms post-entry
                            thisExitSpeedsFullLength(thisExitPreExitLogInd & ~ inClusterLogInd') = NaN;
                            % shorten full length speeds back to the window of interest
                            if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                                thisExitSpeeds(beforeStartFrameNum+1:end) = thisExitSpeedsFullLength(thisExitStartIdx:thisExitEndIdx);
                            elseif exist('afterEndFrameNum','var')
                                thisExitSpeeds(1:(end-afterEndFrameNum)) = thisExitSpeedsFullLength(thisExitStartIdx:thisExitEndIdx);
                            else
                                thisExitSpeeds = thisExitSpeedsFullLength(thisExitStartIdx:thisExitEndIdx);
                            end
                        end
                        % add data to speed matrix
                        exitSpeeds(exitCtr,:) = thisExitSpeeds;
                    end
                    % clear variables
                    clear beforeStartFrameNum
                    clear afterEndFrameNum
                end
            end
            % set maximum speed
            exitSpeeds(abs(exitSpeeds)>1500) = NaN;
            exitSpeeds(exitSpeeds==0) = NaN;
            % smooth speeds
            smoothExitSpeeds = smoothdata(exitSpeeds,2,'movmean',smoothWindow,'includenan');
            % set maximum speed
            smoothExitSpeeds(abs(smoothExitSpeeds)>1500) = NaN;
            
            if ~useManualEvents
                %% save data for pooling
                allEntrySpeeds{fileCtr} = entrySpeeds;
                allExitSpeeds{fileCtr} = exitSpeeds;
                allSmoothEntrySpeeds{fileCtr} = smoothEntrySpeeds;
                allSmoothExitSpeeds{fileCtr} = smoothExitSpeeds;
            end
            
        end
        
        assert(entryCtr == totalEntry+1 & exitCtr == totalExit+1)
        
        if ~useManualEvents
            % pool data across all movies
            allEntrySpeeds = vertcat(allEntrySpeeds{:});
            allExitSpeeds = vertcat(allExitSpeeds{:});
            allSmoothEntrySpeeds = vertcat(allSmoothEntrySpeeds{:});
            allSmoothExitSpeeds = vertcat(allSmoothExitSpeeds{:});
        else
            % rename variables
            allEntrySpeeds = entrySpeeds;
            allExitSpeeds = exitSpeeds;
            allSmoothEntrySpeeds = smoothEntrySpeeds;
            allSmoothExitSpeeds = smoothExitSpeeds;
        end
        
        %% plotting and saving data
        
        % save data
        if useManualEvents
            filename = './figures/entryExitSpeeds/allSmoothEntryExitSpeeds_manualEvent.mat';
        else
            if saveResults
                if enforceInClusterAfterEntryBeforeExit
                    filename = './figures/entryExitSpeeds/allSmoothEntryExitSpeeds_halfFiltered.mat';
                else
                    filename = './figures/entryExitSpeeds/allSmoothEntryExitSpeeds.mat';
                end
            end
        end
        save(filename,'allSmoothEntrySpeeds','allSmoothExitSpeeds','allEntrySpeeds','allExitSpeeds','timeSeries')
        
        if useManualEvents
            
            % entry speed plot depicting single trajectories
            numGraphs = round(totalEntry/maxTrajPerGraph);
            lineColors = distinguishable_colors(totalEntry);
            for graphCtr = 1:numGraphs % go through each graph (each with specified max num of traj plotted)
                entrySpeedsFig = figure; hold on
                set(0, 'CurrentFigure',entrySpeedsFig)
                startTrajIdx = 1+(graphCtr-1)*maxTrajPerGraph;
                if graphCtr*maxTrajPerGraph <= totalEntry
                    endTrajIdx =  graphCtr*maxTrajPerGraph;
                else
                    endTrajIdx = totalEntry;
                end
                for trajCtr = startTrajIdx : endTrajIdx
                    plot(timeSeries.entry,allSmoothEntrySpeeds(trajCtr,:),'color',lineColors(trajCtr,:))
                end
                title('cluster entry speeds')
                xlabel('frames')
                ylabel('speed(microns/s)')
                xlim([timeSeries.entry(1)-20 abs(timeSeries.entry(1)-20)])
                ylim([-500 500])
                legend(entryLegend{startTrajIdx:endTrajIdx})
                figurename = (['figures/entryExitSpeeds/entrySpeedsManualEvents_' strain '_' phase '_graph' num2str(graphCtr)]);
                if saveResults
                    exportfig(entrySpeedsFig,[figurename '.eps'],exportOptions)
                    system(['epstopdf ' figurename '.eps']);
                    system(['rm ' figurename '.eps']);
                end
            end
            
            % exit speed plot depicting single trajectories
            numGraphs = round(totalExit/maxTrajPerGraph);
            lineColors = distinguishable_colors(totalExit);
            for graphCtr = 1:numGraphs % go through each graph (each with specified max num of traj plotted)
                exitSpeedsFig = figure; hold on
                set(0, 'CurrentFigure',exitSpeedsFig)
                startTrajIdx = 1+(graphCtr-1)*maxTrajPerGraph;
                if graphCtr*maxTrajPerGraph <= totalExit
                    endTrajIdx =  graphCtr*maxTrajPerGraph;
                else
                    endTrajIdx = totalExit;
                end
                for trajCtr = startTrajIdx:endTrajIdx
                    plot(timeSeries.exit,allSmoothExitSpeeds(trajCtr,:),'color',lineColors(trajCtr,:))
                end
                title('cluster exit speeds')
                xlabel('frames')
                ylabel('speed(microns/s)')
                xlim([-(timeSeries.exit(end)+20) timeSeries.exit(end)+20])                
                ylim([-500 500])
                legend(exitLegend{startTrajIdx:endTrajIdx},'Location','Northwest')
                figurename = (['figures/entryExitSpeeds/exitSpeedsManualEvents_' strain '_' phase '_graph' num2str(graphCtr)]);
                if saveResults
                    exportfig(exitSpeedsFig,[figurename '.eps'],exportOptions)
                    system(['epstopdf ' figurename '.eps']);
                    system(['rm ' figurename '.eps']);
                end
            end
            
        else
            % mean entry and exit speed plot
            meanEntryExitSpeedsFig = figure; hold on
            set(0,'CurrentFigure',meanEntryExitSpeedsFig)
            plot(timeSeries,nanmean(allSmoothEntrySpeeds,1))
            plot(timeSeries,nanmean(allSmoothExitSpeeds,1))
            title(['mean cluster entry and exit speeds'])
            xlabel('frames')
            ylabel('speed(microns/s)')
            legend('entry','exit')
            if enforceInClusterAfterEntryBeforeExit
                figurename = (['figures/entryExitSpeeds/entryExitSpeedsMeanSmoothed_halfFiltered_' phase]);
            else
                figurename = (['figures/entryExitSpeeds/entryExitSpeedsMeanSmoothed_' phase]);
            end
            if saveResults
                exportfig(meanEntryExitSpeedsFig,[figurename '.eps'],exportOptions)
                system(['epstopdf ' figurename '.eps']);
                system(['rm ' figurename '.eps']);
            end
            
            % mean entry and exit speed plot, separate pos and neg speed
            meanEntryExitSpeeds = figure;
            allSmoothEntrySpeedsPos = allSmoothEntrySpeeds;
            allSmoothEntrySpeedsPos(allSmoothEntrySpeeds<0)=NaN;
            allSmoothEntrySpeedsNeg = allSmoothEntrySpeeds;
            allSmoothEntrySpeedsNeg(allSmoothEntrySpeeds>0)=NaN;
            allSmoothExitSpeedsPos = allSmoothExitSpeeds;
            allSmoothExitSpeedsPos(allSmoothExitSpeeds<0)=NaN;
            allSmoothExitSpeedsNeg = allSmoothExitSpeeds;
            allSmoothExitSpeedsNeg(allSmoothExitSpeeds>0)=NaN;
            subplot(2,1,1); hold on
            plot(timeSeries,nanmean(allSmoothEntrySpeedsPos,1))
            plot(timeSeries,nanmean(allSmoothExitSpeedsPos,1))
            title(['mean positive cluster entry and exit speeds'])
            xlabel('frames')
            ylabel('speed(microns/s)')
            ylim([50 150])
            legend('entry','exit')
            subplot(2,1,2); hold on
            plot(timeSeries,nanmean(allSmoothEntrySpeedsNeg,1))
            plot(timeSeries,nanmean(allSmoothExitSpeedsNeg,1))
            title(['mean negative cluster entry and exit speeds'])
            xlabel('frames')
            ylabel('speed(microns/s)')
            ylim([-150 -50])
            legend('entry','exit')
            if enforceInClusterAfterEntryBeforeExit
                figurename = (['figures/entryExitSpeeds/entryExitSpeedsMeanSmoothedSigned_halfFiltered_' phase]);
            else
                figurename = (['figures/entryExitSpeeds/entryExitSpeedsMeanSmoothedSigned_' phase]);
            end
            if saveResults
                exportfig(meanEntryExitSpeeds,[figurename '.eps'],exportOptions)
                system(['epstopdf ' figurename '.eps']);
                system(['rm ' figurename '.eps']);
            end
            
            % errorbar entry plot
            figure;
            shadedErrorBar(timeSeries,nanmean(allSmoothEntrySpeeds,1),nanstd(allSmoothEntrySpeeds,1),'b');
            title(['mean cluster entry speeds'])
            xlabel('frames')
            ylabel('speed(microns/s)')
            ylim([-150 250])
            if enforceInClusterAfterEntryBeforeExit
                figurename = (['figures/entryExitSpeeds/entrySpeedsMeanErrorSmoothed_halfFiltered_' phase]);
            else
                figurename = (['figures/entryExitSpeeds/entrySpeedsMeanErrorSmoothed_' phase]);
            end
            if saveResults
                entryErrorBarFig = gcf;
                exportfig(entryErrorBarFig,[figurename '.eps'],exportOptions)
                system(['epstopdf ' figurename '.eps']);
                system(['rm ' figurename '.eps']);
            end
            
            % errorbar exit plot
            figure;
            shadedErrorBar(timeSeries,nanmean(allSmoothExitSpeeds,1),nanstd(allSmoothExitSpeeds,1),'r');
            title(['mean cluster exit speeds'])
            xlabel('frames')
            ylabel('speed(microns/s)')
            ylim([-150 250])
            if enforceInClusterAfterEntryBeforeExit
                figurename = (['figures/entryExitSpeeds/exitSpeedsMeanErrorSmoothed_halfFiltered_' phase]);
            else
                figurename = (['figures/entryExitSpeeds/exitSpeedsMeanErrorSmoothed_' phase]);
            end
            if saveResults
                exitErrorBarFig = gcf;
                exportfig(exitErrorBarFig,[figurename '.eps'],exportOptions)
                system(['epstopdf ' figurename '.eps']);
                system(['rm ' figurename '.eps']);
            end
        end
    end
end