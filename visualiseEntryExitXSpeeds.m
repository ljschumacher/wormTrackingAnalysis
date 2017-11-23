clear
close all

%% to improve the script: limit how far to linearly interpolate within naninterp i.e. maybe specify maxinum interpolation window to be 9 or something and leave values NaN if no values nearby;
%% also turn analysis needs debugging

% Script plots midbody speed for worms entering or leaving a cluster using the manually joined red trajectories, with
% the point of entry/exit aligned.

%% set parameters
phase = 'joining'; % 'fullMovie', 'joining', or 'sweeping'.
strains = {'npr1'}; % {'npr1','N2'}
wormnums = {'40'};% {'40'};
preExitDuration = 20; % duration (in seconds) before a worm exits a cluster to be included in the leave cluster analysis
postExitDuration = 20; % duration (in seconds) after a worm exits a cluster to be included in the leave cluster analysis
smoothWindow = 9; % number of frames to smooth over for later trajectory-specific midbody speed calculations
saveResults = true;
alignExitWithStart = true;
addSampleSingleTraj = true;

useManualEvents = true; % manual events only available for joining phase
if useManualEvents
    manualEventMaxDuration = 400; % max number of frames that contains the beginning and end of an annotated event
    maxTrajPerGraph = 10;  
else
    minInOutClusterFrameNum = 5; % unless using manually labelled events, then enter the min number of frames for in/out cluster status to be sustained for an event to be considered
    enforceInClusterAfterEntryBeforeExit = true;
end

midbodyIndcs = 19:33;
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
maxBlobSize_r = 2.5e5;
minSkelLength_r = 850;
maxSkelLength_r = 1500;
minNeighbrDist = 2000;
inClusterNeighbourNum = 3;

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',20,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',15,...
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
                [~,~,annotations] = xlsread(['/data2/shared/data/twoColour/entryExitEvents_' strains{strainCtr} '_' wormnum '_' phase '.xlsx'],1,'A2:I80','basic');
            end
            % count the total number of entry and exit events
            totalEntry = 0;
            totalExit = 0;
            for eventCtr = 1:size(annotations,1)
                if strcmp(annotations{eventCtr,5},'enter')
                    totalEntry = totalEntry + 1;
                elseif strcmp(annotations{eventCtr,5},'exit')
                    totalExit = totalExit+1;
                end
            end
            % pre-allocate matrix to write speeds for each event
            frameRate = 9; % should be the case for most if not all files..
            entrySpeeds = NaN(totalEntry,frameRate*(preExitDuration+postExitDuration)+1+manualEventMaxDuration); % assuming no more than 400 frames for the actual entry/exit event
            exitSpeeds = NaN(totalExit,frameRate*(preExitDuration+postExitDuration)+1+manualEventMaxDuration);
            headAngTotal.entry = NaN(1,totalEntry);
            headAngTotal.exit = NaN(1,totalExit);
            headAngNorm.entry = NaN(1,totalEntry);
            headAngNorm.exit = NaN(1,totalExit);
        end
        
        % initialise counters
        entryCtr = 1;
        exitCtr = 1;
        exitDurations = NaN(1,totalExit); 
        
        %% go through individual movies
        for fileCtr = 1:numFiles
            %% load data
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton'); % in pixels
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            features = h5read(strrep(filename,'skeletons','feat_manual'),'/features_timeseries');
            
            %% generate time series (x), setting time to be zero at the point of entry start (when a worm begins to enter a cluster) or exit end (when a worm fully exits a cluster)
            if useManualEvents
                timeSeries.entry = [-preExitDuration*frameRate:postExitDuration*frameRate+manualEventMaxDuration]/frameRate;
                if alignExitWithStart
                    timeSeries.exit = [-preExitDuration*frameRate:manualEventMaxDuration+postExitDuration*frameRate]/frameRate;
                else
                    timeSeries.exit = [-preExitDuration*frameRate-manualEventMaxDuration:postExitDuration*frameRate]/frameRate;
                end

            else
                timeSeries = [-preExitDuration*frameRate:postExitDuration*frameRate]/frameRate;
            end
            
            if ~useManualEvents
                
                %% calculate midbody signed speed (from reversalAnalysisBodyWall.m)
                [midbodySpeed,~,~,midbodySpeedSigned] = calculateSpeedsFromSkeleton(trajData,skelData,...
                    midbodyIndcs,pixelsize,frameRate,false,smoothWindow);
                % ignore first and last frames of each worm's track
                wormChangeIndcs = gradient(double(trajData.worm_index_manual))~=0;
                midbodySpeedSigned(wormChangeIndcs)=NaN;
            else
                % load unsorted xy coords; sort later using worm index for interpolation
                xcoords = squeeze(skelData(1,:,:)); 
                ycoords = squeeze(skelData(2,:,:));
            end
            
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
            trajData.filtered = trajData.filtered...
                & filterSkelLength(skelData,pixelsize,minSkelLength_r,maxSkelLength_r);
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
            
            %% entry case with manually labelled events
            if useManualEvents
                % go through each annotated event
                for eventCtr = 1:size(annotations,1)
                    
                    %% find entry events with matching recording file name to the file currently under consideration
                    recordingNumber = annotations{eventCtr,1};
                    recordingNumber = recordingNumber(2:end); % remove the 'r' before the number
                    
                    if contains(filename, recordingNumber) & strcmp(annotations{eventCtr,5},'enter')
                        %% collect information on worm index and entry frames
                        wormIndex = annotations{eventCtr,2};
                        % get the annotated entry start frame
                        thisEntryStartFrame = annotations{eventCtr,7}; 
                         % get the annotated entry finish frame
                        thisEntryEndFrame= annotations{eventCtr,8};
                        % trim hand annotation in case event end goes beyond phase end
                        if thisEntryEndFrame > lastPhaseFrame 
                            thisEntryEndFrame = lastPhaseFrame;
                        end
                        if thisEntryEndFrame-thisEntryStartFrame +1 > manualEventMaxDuration
                            warning('event duration exceeds allocated maximum duration, increase value for manualEventMaxDuration variable')
                        end
                        
                        %% extend for a specified duration before and after the entry point
                        thisEntryXStartFrame = thisEntryStartFrame - preExitDuration*frameRate;
                        % check that the extended start doesn't go beyond the phase of interest
                        if thisEntryXStartFrame <= firstPhaseFrame
                            % take note of omitted frames for alignment purposes
                            beforeStartFrameNum = firstPhaseFrame-thisEntryXStartFrame; 
                             % exclude frames before the start of the specified phase
                            thisEntryXStartFrame = firstPhaseFrame;
                        end
                        thisEntryXEndFrame = thisEntryEndFrame + postExitDuration*frameRate;
                         % check that the extended start doesn't go beyond the phase of interest
                        if  thisEntryXEndFrame > lastPhaseFrame
                            % take note of omitted frames for alignment purposes
                            afterEndFrameNum =  thisEntryXEndFrame -lastPhaseFrame; 
                            % exclude frames beyond the end of the specified phase
                            thisEntryXEndFrame = lastPhaseFrame; 
                        end
                        
                        %% get aligned list of frames for the event
                        thisEntrySpeeds = NaN(1,frameRate*(preExitDuration+postExitDuration)+1+manualEventMaxDuration);
                        thisEntryXFrames = NaN(1,frameRate*(preExitDuration+postExitDuration)+1+manualEventMaxDuration);
                        thisEntryHeadSpeedLogInd = false(1,frameRate*(preExitDuration+postExitDuration)+1+manualEventMaxDuration);
                        if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                            thisEntryXFrames(beforeStartFrameNum+1:end) = thisEntryXStartFrame:thisEntryXEndFrame;
                        elseif exist('afterEndFrameNum','var')
                            thisEntryXFrames(1:end-afterEndFrameNum) = thisEntryXStartFrame:thisEntryXEndFrame;
                        else
                            thisEntryXFrames = thisEntryXStartFrame:thisEntryXEndFrame;
                        end

                        %% sort midbody speed and frames for the worm at consideration
                        % get the indices for the worm of interest
                        uniqueWormIdx = find(trajData.worm_index_manual == wormIndex); 
                        % extract skeleton data for worm of interest
                        xcoordsSorted = xcoords(:,uniqueWormIdx); 
                        ycoordsSorted = ycoords(:,uniqueWormIdx);
                        frameNumberSorted = trajData.frame_number(uniqueWormIdx);
                        
                        %% interpolate over NaN values for sorted xy coordinates
                        for nodeCtr = 1:size(xcoords,1)
                            xcoordsNode = xcoordsSorted(nodeCtr,:);
                            ycoordsNode = ycoordsSorted(nodeCtr,:);
                            xcoordsNode = naninterp(xcoordsNode); % naninterp only works for vectors so go node by node
                            ycoordsNode = naninterp(ycoordsNode);
                            xcoordsSorted(nodeCtr,:) = xcoordsNode;
                            ycoordsSorted(nodeCtr,:) = ycoordsNode;
                        end
                        
                        %% calculate midbodyspeed using sorted, interpolated xy coordinates
                        % centroids of midbody skeleton
                        x = mean(xcoordsSorted(midbodyIndcs,:))*pixelsize;
                        y = mean(ycoordsSorted(midbodyIndcs,:))*pixelsize;
                        % change in centroid position over time
                        dxdt = gradient(x)*frameRate;
                        dydt = gradient(y)*frameRate;
                        % speed and velocity
                        dFramedt = gradient(double(frameNumberSorted))';
                        midbodySpeed = sqrt(dxdt.^2 + dydt.^2)./dFramedt;
                        velocity_x = dxdt./dFramedt;
                        velocity_y = dydt./dFramedt;
                        % signed speed calculation
                        % direction of segments pointing along midbody
                        [~, dxds] = gradient(xcoordsSorted(midbodyIndcs,:),-1);
                        [~, dyds] = gradient(ycoordsSorted(midbodyIndcs,:),-1);
                        % sign speed based on relative orientation of velocity to body
                        midbodySpeedSigned = getSignedSpeed([velocity_x; velocity_y],[mean(dxds); mean(dyds)]);

                        %% go through each frame to obtain speed
                        for frameCtr = 1:length(thisEntryXFrames)
                            frameNumber = thisEntryXFrames(frameCtr);
                            if ~isnan(frameNumber)
                                wormFrameLogInd = frameNumberSorted == frameNumber;
                                if nnz(wormFrameLogInd)~=0
                                    assert(nnz(wormFrameLogInd) ==1);
                                    % obtain speed
                                    thisEntrySpeeds(frameCtr) = midbodySpeedSigned(wormFrameLogInd);
                                    % obtain xy coordinates for subsequent head angle calculation
                                    if frameNumber<max(thisEntryXFrames)-postExitDuration*frameRate % only interested in the time points before entry point
                                        xcoordsSortedWorm(:,frameCtr) = xcoordsSorted(:,wormFrameLogInd);
                                        ycoordsSortedWorm(:,frameCtr) = ycoordsSorted(:,wormFrameLogInd);
                                    end
                                end
                            end
                        end
                        % remove 0 values and rename xy coordinate variable
                        if exist('xcoordsSortedWorm','var')
                            xcoordsSortedWorm(xcoordsSortedWorm==0)=[];
                            ycoordsSortedWorm(ycoordsSortedWorm==0)=[];
                            xcoordsSorted = reshape(xcoordsSortedWorm,49,[])';
                            ycoordsSorted = reshape(ycoordsSortedWorm,49,[])';
                            clear xcoordsSortedWorm
                            clear ycoordsSortedWorm
                        else
                            warning(['ignore data for entryCtr = ' num2str(entryCtr)])
                        end
                        
                        %% write speeds
                        entrySpeeds(entryCtr,1:length(thisEntrySpeeds)) = thisEntrySpeeds;
                        % clear variables
                        clear beforeStartFrameNum
                        clear afterEndFrameNum
                        % write event number to legend
                        entryLegend{entryCtr} = num2str(eventCtr);
                        
                        %% calculate head angles
                        % calculate head angle changes per frame
                        [headAngleDiff, framesElapsed] = getHeadAngleDiff(xcoordsSorted,ycoordsSorted, 'bodywall', true, frameRate);
                        % calculate total head turn over full trajectory
                        headAngTotal.entry(entryCtr) = abs(nansum(headAngleDiff));
                        % normalise head turn over path length
                        xcoordsSorted = nanmean(xcoordsSorted,2);
                        ycoordsSorted = nanmean(ycoordsSorted,2);
                        xdiff = xcoordsSorted(2:end) - xcoordsSorted(1:end-1);
                        ydiff = ycoordsSorted(2:end) - ycoordsSorted(1:end-1);
                        pathLength = sum(sqrt(xdiff.^2+ydiff.^2)); % calculate path length
                        headAngNorm.entry(entryCtr) = headAngTotal.entry(entryCtr)/pathLength;
                        % calculate head angular speed
                        % headAngSpeed.entry(entryCtr) = headAngTotal.entry(entryCtr)/framesElapsed*frameRate;
                        headAngSpeed.entry{entryCtr} = headAngleDiff*frameRate;
                        % save xy coordinates
                        sortedxcoords.entry{entryCtr} = xcoordsSorted;
                        sortedycoords.entry{entryCtr} = ycoordsSorted;
                        
                        %% update counter
                        entryCtr = entryCtr +1;
                        assert(entryCtr<=totalEntry+1);
                        
                    end
                end
                
                %% without manual event annotation, we need to identify our own events
            else
                
                %% get indices for cluster entry points
                enterClusterPoint = find(enterClusterStartLogInd);
                entrySpeeds = NaN(numel(enterClusterPoint),preExitDuration*frameRate +1+ postExitDuration*frameRate);
                % loop through each entry
                for entryCtr = 1:numel(enterClusterPoint)
                    thisEntryIdx = enterClusterPoint(entryCtr);
                    
                    %% get the worm index
                    wormIndex = trajData.worm_index_manual(thisEntryIdx);
                    
                    %% check that the same worm stays in cluster for minumum number of frames
                    if thisEntryIdx+minInOutClusterFrameNum <= length(trajData.frame_number) &&...
                            nnz(inClusterLogInd(thisEntryIdx:(thisEntryIdx+minInOutClusterFrameNum))) == 1+minInOutClusterFrameNum...
                            && nnz(trajData.worm_index_manual(thisEntryIdx:(thisEntryIdx+minInOutClusterFrameNum)) == wormIndex) == 1+minInOutClusterFrameNum
                        
                        %% expand for a specified number of entries before and after the entry point
                        thisEntryStartIdx = thisEntryIdx-preExitDuration*frameRate;
                        if thisEntryStartIdx <= 0
                            % take note of omitted frames for alignment purposes
                            beforeStartFrameNum = 1-thisEntryStartIdx; 
                            % exclude entries below index 0
                            thisEntryStartIdx = 1; 
                        end
                        thisEntryEndIdx = thisEntryIdx + postExitDuration*frameRate;
                        if thisEntryEndIdx > length(trajData.frame_number)
                             % take note of omitted frames for alignment purposes
                            afterEndFrameNum = thisEntryEndIdx - length(trajData.frame_number);
                            % exclude entries above highest index
                            thisEntryEndIdx = length(trajData.frame_number); 
                        end
                        
                        %% generate logical index for this entry
                        thisEntryLogInd = false(1,length(midbodySpeedSigned));
                        thisEntryLogInd(thisEntryStartIdx:thisEntryEndIdx) = true;
                        
                        %% get midbody signed speeds (y) for this entry
                        thisEntryXSpeeds = NaN(1,size(entrySpeeds,2));
                        % this keeps the alignment of the entries in case they go below or above min/max index
                        if exist('beforeStartFrameNum','var') 
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
                        %% add data to speed matrix
                        entrySpeeds(entryCtr,:) = thisEntryXSpeeds;
                    end
                    % clear variables
                    clear beforeStartFrameNum
                    clear afterEndFrameNum
                end
            end
            
            % set maximum speed and remove 0 speed
            entrySpeeds(abs(entrySpeeds)>1500) = NaN;
            % smooth speeds
            if smoothWindow==0
                smoothEntrySpeeds = entrySpeeds;
            elseif smoothWindow >3
                smoothEntrySpeeds = smoothdata(entrySpeeds,2,'movmean',smoothWindow,'includenan');
            else
                smoothEntrySpeeds = smoothdata(entrySpeeds,2,'movmean',smoothWindow);
            end
            % set maximum speed
            smoothEntrySpeeds(abs(smoothEntrySpeeds)>1500) = NaN;
            
            %% obtain speed over time course for cluster exit events
            
            %% exit cases with manually labelled events
            if useManualEvents
                
                % go through each annotated event
                for eventCtr =1:size(annotations,1)
                    
                    %% find exit events with matching recording file name to the file currently under consideration
                    recordingNumber = annotations{eventCtr,1};
                    recordingNumber = recordingNumber(2:end); % remove the 'r' before the number
                    
                    if contains(filename, recordingNumber) & strcmp(annotations{eventCtr,5},'exit')
                        
                        %% collect information on worm index and exit frames
                        wormIndex = annotations{eventCtr,2};
                        thisExitStartFrame = annotations{eventCtr,7}; % get the annotated entry start frame
                        thisExitEndFrame= annotations{eventCtr,8}; % get the annotated entry finish frame
                        exitDurations(exitCtr) = (thisExitEndFrame - thisExitStartFrame)/frameRate; 
                        if thisExitEndFrame > lastPhaseFrame % trim hand annotation in case event end goes beyond phase end
                            thisExitEndFrame = lastPhaseFrame;
                        end
                        if thisExitEndFrame-thisExitStartFrame +1 > manualEventMaxDuration
                            warning('event duration exceeds allocated maximum duration, increase value for manualEventMaxDuration variable')
                        end
                        
                        %% extend for a specified duration before and after the exit point
                        thisExitXStartFrame = thisExitStartFrame - preExitDuration*frameRate;
                        % check that the extended start doesn't go beyond the phase of interest
                        if thisExitXStartFrame <= firstPhaseFrame
                            beforeStartFrameNum = firstPhaseFrame-thisExitXStartFrame; % take note of omitted frames for alignment purposes
                            thisExitXStartFrame = firstPhaseFrame; % exclude frames before the start of the specified phase
                        end
                        thisExitXEndFrame = thisExitEndFrame + postExitDuration*frameRate;
                        % check that the extended start doesn't go beyond the phase of interest
                        if  thisExitXEndFrame > lastPhaseFrame
                            afterEndFrameNum =  thisExitXEndFrame -lastPhaseFrame; % take note of omitted frames for alignment purposes
                            thisExitXEndFrame = lastPhaseFrame; % exclude frames beyond the end of the specified phase
                        end
                        
                        %% get aligned list of frames for the event (align speed vectors using the end of the exit event i.e. when a worm fully exits a cluster)
                        exitNumFrames = thisExitEndFrame - thisExitStartFrame;
                        exitDurationFiller = manualEventMaxDuration - exitNumFrames ; % number of empty frames to add to the start/end of speed vector to keep alignment for end/start of exit
                        thisExitSpeeds = NaN(1,frameRate*(preExitDuration+postExitDuration)+1+manualEventMaxDuration);
                        thisExitXFrames = NaN(1,frameRate*(preExitDuration+postExitDuration)+1+manualEventMaxDuration);
                        
                        if alignExitWithStart
                            if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                                thisExitXFrames(beforeStartFrameNum+1:end-exitDurationFiller) = thisExitXStartFrame:thisExitXEndFrame;
                            elseif exist('afterEndFrameNum','var')
                                thisExitXFrames(1:(end-afterEndFrameNum-exitDurationFiller)) = thisExitXStartFrame:thisExitXEndFrame;
                            else
                                thisExitXFrames(1:end-exitDurationFiller) = thisExitXStartFrame:thisExitXEndFrame;
                            end
                        else
                            if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                                thisExitXFrames(beforeStartFrameNum+exitDurationFiller+1:end) = thisExitXStartFrame:thisExitXEndFrame;
                            elseif exist('afterEndFrameNum','var')
                                thisExitXFrames(exitDurationFiller+1:(end-afterEndFrameNum)) = thisExitXStartFrame:thisExitXEndFrame;
                            else
                                thisExitXFrames(exitDurationFiller+1:end) = thisExitXStartFrame:thisExitXEndFrame;
                            end
                        end
                        
                        %% sort midbody speed and frames for the worm at consideration
                        % get the indices for the worm of interest
                        uniqueWormIdx = find(trajData.worm_index_manual == wormIndex);
                        % extract skeleton data for worm of interest
                        xcoordsSorted = xcoords(:,uniqueWormIdx);
                        ycoordsSorted = ycoords(:,uniqueWormIdx);
                        frameNumberSorted = trajData.frame_number(uniqueWormIdx);
                        
                        %% interpolate over NaN values for sorted xy coordinates
                        for nodeCtr = 1:size(xcoords,1)
                            xcoordsNode = xcoordsSorted(nodeCtr,:);
                            ycoordsNode = ycoordsSorted(nodeCtr,:);
                            xcoordsNode = naninterp(xcoordsNode);%naninterp only works for vectors so go node by node
                            ycoordsNode = naninterp(ycoordsNode);
                            xcoordsSorted(nodeCtr,:) = xcoordsNode;
                            ycoordsSorted(nodeCtr,:) = ycoordsNode;
                        end
                        
                        %% calculate midbodyspeed using sorted, interpolated xy coordinates
                        % centroids of midbody skeleton
                        x = mean(xcoordsSorted(midbodyIndcs,:))*pixelsize;
                        y = mean(ycoordsSorted(midbodyIndcs,:))*pixelsize;
                        % change in centroid position over time
                        dxdt = gradient(x)*frameRate;
                        dydt = gradient(y)*frameRate;
                        % speed and velocity
                        dFramedt = gradient(double(frameNumberSorted))';
                        midbodySpeed = sqrt(dxdt.^2 + dydt.^2)./dFramedt;
                        velocity_x = dxdt./dFramedt;
                        velocity_y = dydt./dFramedt;
                        % signed speed calculation
                        % direction of segments pointing along midbody
                        [~, dxds] = gradient(xcoordsSorted,-1);
                        [~, dyds] = gradient(ycoordsSorted,-1);
                        % sign speed based on relative orientation of velocity to body
                        midbodySpeedSigned = getSignedSpeed([velocity_x; velocity_y],[mean(dxds); mean(dyds)]);
                        
                        %% go through each frame
                        for frameCtr = 1:length(thisExitXFrames)
                            frameNumber = thisExitXFrames(frameCtr);
                            if ~isnan(frameNumber)
                                wormFrameLogInd = frameNumberSorted == frameNumber;
                                if nnz(wormFrameLogInd)~=0
                                    assert(nnz(wormFrameLogInd) ==1);
                                    % obtain speed
                                    thisExitSpeeds(frameCtr) = midbodySpeedSigned(wormFrameLogInd);
                                    % obtain xy coordinates for subsequent head angle calculation
                                    if frameNumber>max(thisExitXFrames)-postExitDuration*frameRate-1 % only interested in coordinates for the time points after exit point
                                        xcoordsSortedWorm(:,frameCtr) = xcoordsSorted(:,wormFrameLogInd);
                                        ycoordsSortedWorm(:,frameCtr) = ycoordsSorted(:,wormFrameLogInd);
                                    end
                                end
                            end
                        end
                        
                        % remove 0 values and rename xy coordinate variable
                        if exist('xcoordsSortedWorm','var')
                        xcoordsSortedWorm(xcoordsSortedWorm==0)=[];
                        ycoordsSortedWorm(ycoordsSortedWorm==0)=[];
                        xcoordsSorted = reshape(xcoordsSortedWorm,49,[])';
                        ycoordsSorted = reshape(ycoordsSortedWorm,49,[])';
                        clear xcoordsSortedWorm
                        clear ycoordsSortedWorm
                        else
                            warning(['ignore data for exitCtr = ' num2str(exitCtr)])
                        end
                        
                        %% write speeds
                        exitSpeeds(exitCtr,1:length(thisExitSpeeds)) = thisExitSpeeds;
                        % clear variables
                        clear beforeStartFrameNum
                        clear afterEndFrameNum
                        % write event number to legend
                        exitLegend{exitCtr} = num2str(eventCtr);
                        
%                         %% calculate head angles
% 
%                         % calculate head angle changes per frame
%                         [headAngleDiff, framesElapsed] = getHeadAngleDiff(xcoordsSorted,ycoordsSorted, 'bodywall', true, frameRate);
%                         % calculate total head turn over full trajectory
%                         headAngTotal.exit(exitCtr) = abs(nansum(headAngleDiff));
%                         % normalise head turn over path length
%                         xcoordsSorted = nanmean(xcoordsSorted,2);
%                         ycoordsSorted = nanmean(ycoordsSorted,2);
%                         xdiff = xcoordsSorted(2:end) - xcoordsSorted(1:end-1);
%                         ydiff = ycoordsSorted(2:end) - ycoordsSorted(1:end-1);
%                         pathLength = sum(sqrt(xdiff.^2+ydiff.^2)); % calculate path length
%                         headAngNorm.exit(exitCtr) = headAngTotal.exit(exitCtr)/pathLength;
%                         % calculate head angular speed
%                         %headAngSpeed.exit(exitCtr) = headAngTotal.exit(exitCtr)/framesElapsed*frameRate;
%                         headAngSpeed.exit{exitCtr} = headAngleDiff*frameRate;
%                         % save xy coordinates
%                         sortedxcoords.exit{exitCtr} = xcoordsSorted;
%                         sortedycoords.exit{exitCtr} = ycoordsSorted;
                        
                        %% update exit counter
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
            % smooth speeds
            if smoothWindow ==0
                smoothExitSpeeds = exitSpeeds;
            elseif smoothWindow >3
                smoothExitSpeeds = smoothdata(exitSpeeds,2,'movmean',smoothWindow,'includenan');
            else
                smoothExitSpeeds = smoothdata(exitSpeeds,2,'movmean',smoothWindow);
            end
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
        
        
        if useManualEvents
            % rename variables (so they are the same as the ~useManualEvents case)
            allEntrySpeeds = entrySpeeds;
            allExitSpeeds = exitSpeeds;
            allSmoothEntrySpeeds = smoothEntrySpeeds;
            allSmoothExitSpeeds = smoothExitSpeeds;
        else
            % pool data across all movies
            allEntrySpeeds = vertcat(allEntrySpeeds{:});
            allExitSpeeds = vertcat(allExitSpeeds{:});
        end
        
        %% plotting and saving data
        
%         % save data
%         if saveResults
%             if useManualEvents
%                 filename = './figures/entryExitSpeeds/allSmoothEntryExitSpeeds_manualEvent.mat';
%                 save(filename,'allSmoothEntrySpeeds','allSmoothExitSpeeds','allEntrySpeeds','allExitSpeeds','timeSeries',...
%                     'headAngTotal','headAngNorm','headAngSpeed','sortedxcoords','sortedycoords')
%             else
%                 if enforceInClusterAfterEntryBeforeExit
%                     filename = './figures/entryExitSpeeds/allSmoothEntryExitSpeeds_halfFiltered.mat';
%                 else
%                     filename = './figures/entryExitSpeeds/allSmoothEntryExitSpeeds.mat';
%                 end
%                 save(filename,'allSmoothEntrySpeeds','allSmoothExitSpeeds','allEntrySpeeds','allExitSpeeds','timeSeries')
%             end
%         end
        
        
        if useManualEvents
            
%             % entry speed plot depicting single trajectories
%             numGraphs = round(totalEntry/maxTrajPerGraph);
%             lineColors = distinguishable_colors(totalEntry);
%             for graphCtr = 1:numGraphs % go through each graph (each with specified max num of traj plotted)
%                 entrySpeedsFig = figure; hold on
%                 set(0, 'CurrentFigure',entrySpeedsFig)
%                 startTrajIdx = 1+(graphCtr-1)*maxTrajPerGraph;
%                 if graphCtr*maxTrajPerGraph <= totalEntry
%                     endTrajIdx =  graphCtr*maxTrajPerGraph;
%                 else
%                     endTrajIdx = totalEntry;
%                 end
%                 for trajCtr = startTrajIdx : endTrajIdx
%                     plot(timeSeries.entry,allSmoothEntrySpeeds(trajCtr,:),'color',lineColors(trajCtr,:))
%                 end
%                 vline(0,'k')
%                 title('cluster entry speeds')
%                 xlabel('time(s)')
%                 ylabel('speed(microns/s)')
%                 xlim([-preExitDuration postExitDuration])
%                 ylim([-500 500])
%                 legend(entryLegend{startTrajIdx:endTrajIdx})
%                 figurename = (['figures/entryExitSpeeds/entrySpeedsManualEvents_' strain '_' phase '_graph' num2str(graphCtr)  '_smoothWindow' num2str(smoothWindow)]);
%                 if saveResults
%                     exportfig(entrySpeedsFig,[figurename '.eps'],exportOptions)
%                     system(['epstopdf ' figurename '.eps']);
%                     system(['rm ' figurename '.eps']);
%                 end
%             end
%             
%             % exit speed plot depicting single trajectories
%             numGraphs = round(totalExit/maxTrajPerGraph);
%             lineColors = distinguishable_colors(totalExit);
%             for graphCtr = 1:numGraphs % go through each graph (each with specified max num of traj plotted)
%                 exitSpeedsFig = figure; hold on
%                 set(0, 'CurrentFigure',exitSpeedsFig)
%                 startTrajIdx = 1+(graphCtr-1)*maxTrajPerGraph;
%                 if graphCtr*maxTrajPerGraph <= totalExit
%                     endTrajIdx =  graphCtr*maxTrajPerGraph;
%                 else
%                     endTrajIdx = totalExit;
%                 end
%                 for trajCtr = startTrajIdx:endTrajIdx
%                     plot(timeSeries.exit,allSmoothExitSpeeds(trajCtr,:),'color',lineColors(trajCtr,:))
%                 end
%                 vline(0,'k')
%                 title('cluster exit speeds')
%                 xlabel('time(s)')
%                 ylabel('speed(microns/s)')
%                 xlim([-preExitDuration postExitDuration])
%                 ylim([-500 500])
%                 legend(exitLegend{startTrajIdx:endTrajIdx},'Location','Northwest')
%                 figurename = (['figures/entryExitSpeeds/exitSpeedsManualEvents_' strain '_' phase '_graph' num2str(graphCtr) '_smoothWindow' num2str(smoothWindow)]);
%                 if saveResults
%                     exportfig(exitSpeedsFig,[figurename '.eps'],exportOptions)
%                     system(['epstopdf ' figurename '.eps']);
%                     system(['rm ' figurename '.eps']);
%                 end
%             end
%             
%             % average entry and exit speeds plot (of absolute value)
%             figure;
%             % take absolute value
%             allSmoothEntrySpeeds = abs(allSmoothEntrySpeeds);
%             allSmoothExitSpeeds = abs(allSmoothExitSpeeds);
%             % entry subplot
%             subplot(1,2,1); hold on
%             shadedErrorBar(timeSeries.entry,nanmean(allSmoothEntrySpeeds,1),nanstd(allSmoothEntrySpeeds,1),'c');
%             title('cluster entry')
%             xlim([-preExitDuration postExitDuration])
%             xticks(linspace(-preExitDuration,postExitDuration,9))
%             ylim([0 350])
%             xlabel('time(s)')
%             ylabel('speed (microns/s)')
%             legend(['n=' num2str(totalEntry)])
%             vline(0,'k')
%             % exit subplot
%             subplot(1,2,2); hold on
%             shadedErrorBar(timeSeries.exit,nanmean(allSmoothExitSpeeds,1),nanstd(allSmoothExitSpeeds,1),'m');
%             title('cluster exit')
%             xlim([-preExitDuration postExitDuration])
%             xticks(linspace(-preExitDuration,postExitDuration,9))
%             ylim([0 350])
%             xlabel('time(s)')
%             ylabel('speed (\mum/s)')
%             legend(['n=' num2str(totalExit)])
%             vline(0,'k')
%             if ~alignExitWithStart
%                 exitDuration = median(exitDurations);
%                 vline(-exitDuration,'g') % draw a reference exit start line based on median exit duration before worms fully exiting cluster
%             end
%             %
%             figurename = (['figures/entryExitSpeeds/entryExitSpeedsManualEvents_' strain '_' phase '_avgGraph_smoothWindow' num2str(smoothWindow)]);
%             if alignExitWithStart
%                 figurename = ([figurename '_alignWithStart']);
%             else
%                 figurename = ([figurename '_alignWithEnd']);
%             end
%             if saveResults
%                 avgEntryExitSpeedFig = gcf;
%                 exportfig(avgEntryExitSpeedFig,[figurename '.eps'],exportOptions)
%                 system(['epstopdf ' figurename '.eps']);
%                 system(['rm ' figurename '.eps']);
%             end
            
            % entry stand-alone plot
            figure; hold on
            shadedErrorBar(timeSeries.entry,nanmean(allSmoothEntrySpeeds,1),nanstd(allSmoothEntrySpeeds,1),'k');
            if addSampleSingleTraj
                load('/data2/shared/data/twoColour/event49.mat')
                plot(timeSeries.entry,thisEventSpeeds,'r');
            end
            title('cluster entry')
            xlim([-preExitDuration postExitDuration])
            xticks(linspace(-preExitDuration,postExitDuration,9))
            ylim([0 350])
            xlabel('time(s)')
            ylabel('speed (\mum/s)')
            legend(['n=' num2str(totalEntry)])
            vline(0,'k')
            if addSampleSingleTraj
                figurename = (['figures/entryExitSpeeds/entrySpeedsManualEvents_' strain '_' phase '_avgEntryGraph_plusEvent49_smoothWindow' num2str(smoothWindow)]);
            else
                figurename = (['figures/entryExitSpeeds/entrySpeedsManualEvents_' strain '_' phase '_avgEntryGraph_smoothWindow' num2str(smoothWindow)]);
            end
            if saveResults
                avgEntrySpeedFig = gcf;
                exportfig(avgEntrySpeedFig,[figurename '.eps'],exportOptions)
                system(['epstopdf ' figurename '.eps']);
                system(['rm ' figurename '.eps']);
            end
%             
%             % exit stand-alone plot
%             figure;
%             shadedErrorBar(timeSeries.exit,nanmean(allSmoothExitSpeeds,1),nanstd(allSmoothExitSpeeds,1),'k');
%             title('cluster exit')
%             xlim([-preExitDuration postExitDuration])
%             xticks(linspace(-preExitDuration,postExitDuration,9))
%             ylim([0 350])
%             xlabel('time(s)')
%             ylabel('speed (\mum/s)')
%             legend(['n=' num2str(totalExit)])
%             vline(0,'k')
%             if ~alignExitWithStart
%                 exitDuration = median(exitDurations);
%                 vline(-exitDuration,'g') % draw a reference exit start line based on median exit duration before worms fully exiting cluster
%             end
%             figurename = (['figures/entryExitSpeeds/exitSpeedsManualEvents_' strain '_' phase '_avgExitGraph_smoothWindow' num2str(smoothWindow)]);
%             if alignExitWithStart
%                 figurename = ([figurename '_alignWithStart']);
%             else
%                 figurename = ([figurename '_alignWithEnd']);
%             end
%             if saveResults
%                 avgExitSpeedFig = gcf;
%                 exportfig(avgExitSpeedFig,[figurename '.eps'],exportOptions)
%                 system(['epstopdf ' figurename '.eps']);
%                 system(['rm ' figurename '.eps']);
%             end
%             
%             % entry and exit head angular speed plots
%             headAngFig = figure; 
%             %
%             subplot(2,2,1);
%             histogram(headAngTotal.entry,'BinWidth',pi/4,'Normalization','count')
%             title([strains{strainCtr} ' entry'],'FontWeight','normal')
%             xlabel('Total head angle change (radian)')
%             ylabel('Count')
%             xlim([0 7])
%             ylim([0 10])
%             %
%             subplot(2,2,2);
%             headAngSpeedEntry = vertcat(headAngSpeed.entry{:});
%             histogram(headAngSpeedEntry,'BinWidth',pi/4,'Normalization','count')
%             title([strains{strainCtr} ' entry'],'FontWeight','normal')
%             xlabel('Head angular speed (radian/s)')
%             ylabel('Count')
%             xlim([0 7])
%             ylim([0 2500])
%             %
%             subplot(2,2,3);
%             histogram(headAngTotal.exit,'BinWidth',pi/4,'Normalization','count')
%             title([strains{strainCtr} ' exit'],'FontWeight','normal')
%             xlabel('Total head angle change (radian)')
%             ylabel('Count')
%             xlim([0 7])
%             ylim([0 10])
%             %
%             subplot(2,2,4);
%             headAngSpeedExit = vertcat(headAngSpeed.exit{:});
%             histogram(headAngSpeedExit,'BinWidth',pi/4,'Normalization','count')
%             title([strains{strainCtr} ' exit'],'FontWeight','normal')
%             xlabel('Head angular speed (radian/s)')
%             ylabel('Count')
%             xlim([0 7])
%             ylim([0 2500])
%             %
%             set(headAngFig,'PaperUnits','centimeters')
%             figurename = (['figures/entryExitSpeeds/exitSpeedsManualEvents_' strain '_' phase ' headAngles']);
%             if saveResults
%                 exportfig(headAngFig,[figurename '.eps'],exportOptions)
%                 system(['epstopdf ' figurename '.eps']);
%                 system(['rm ' figurename '.eps']);
%             end
%             
%             % not using manually annotated traj
%         else
%             % mean entry and exit speed plot
%             meanEntryExitSpeedsFig = figure; hold on
%             set(0,'CurrentFigure',meanEntryExitSpeedsFig)
%             plot(timeSeries,nanmean(allSmoothEntrySpeeds,1))
%             plot(timeSeries,nanmean(allSmoothExitSpeeds,1))
%             title(['mean cluster entry and exit speeds'])
%             xlabel('frames')
%             ylabel('speed(microns/s)')
%             legend('entry','exit')
%             if enforceInClusterAfterEntryBeforeExit
%                 figurename = (['figures/entryExitSpeeds/entryExitSpeedsMeanSmoothed_halfFiltered_' phase]);
%             else
%                 figurename = (['figures/entryExitSpeeds/entryExitSpeedsMeanSmoothed_' phase]);
%             end
%             if saveResults
%                 exportfig(meanEntryExitSpeedsFig,[figurename '.eps'],exportOptions)
%                 system(['epstopdf ' figurename '.eps']);
%                 system(['rm ' figurename '.eps']);
%             end
%             
%             % mean entry and exit speed plot, separate pos and neg speed
%             meanEntryExitSpeeds = figure;
%             allSmoothEntrySpeedsPos = allSmoothEntrySpeeds;
%             allSmoothEntrySpeedsPos(allSmoothEntrySpeeds<0)=NaN;
%             allSmoothEntrySpeedsNeg = allSmoothEntrySpeeds;
%             allSmoothEntrySpeedsNeg(allSmoothEntrySpeeds>0)=NaN;
%             allSmoothExitSpeedsPos = allSmoothExitSpeeds;
%             allSmoothExitSpeedsPos(allSmoothExitSpeeds<0)=NaN;
%             allSmoothExitSpeedsNeg = allSmoothExitSpeeds;
%             allSmoothExitSpeedsNeg(allSmoothExitSpeeds>0)=NaN;
%             subplot(2,1,1); hold on
%             plot(timeSeries,nanmean(allSmoothEntrySpeedsPos,1))
%             plot(timeSeries,nanmean(allSmoothExitSpeedsPos,1))
%             title(['mean positive cluster entry and exit speeds'])
%             xlabel('frames')
%             ylabel('speed(microns/s)')
%             ylim([50 150])
%             legend('entry','exit')
%             subplot(2,1,2); hold on
%             plot(timeSeries,nanmean(allSmoothEntrySpeedsNeg,1))
%             plot(timeSeries,nanmean(allSmoothExitSpeedsNeg,1))
%             title(['mean negative cluster entry and exit speeds'])
%             xlabel('frames')
%             ylabel('speed(microns/s)')
%             ylim([-150 -50])
%             legend('entry','exit')
%             if enforceInClusterAfterEntryBeforeExit
%                 figurename = (['figures/entryExitSpeeds/entryExitSpeedsMeanSmoothedSigned_halfFiltered_' phase]);
%             else
%                 figurename = (['figures/entryExitSpeeds/entryExitSpeedsMeanSmoothedSigned_' phase]);
%             end
%             if saveResults
%                 exportfig(meanEntryExitSpeeds,[figurename '.eps'],exportOptions)
%                 system(['epstopdf ' figurename '.eps']);
%                 system(['rm ' figurename '.eps']);
%             end
%             
%             % errorbar entry plot
%             
%             figure;
%             shadedErrorBar(timeSeries,nanmean(allSmoothEntrySpeeds,1),nanstd(allSmoothEntrySpeeds,1),'b');
%             title(['mean cluster entry speeds'])
%             xlabel('frames')
%             ylabel('speed(microns/s)')
%             ylim([-150 250])
%             if enforceInClusterAfterEntryBeforeExit
%                 figurename = (['figures/entryExitSpeeds/entrySpeedsMeanErrorSmoothed_halfFiltered_' phase]);
%             else
%                 figurename = (['figures/entryExitSpeeds/entrySpeedsMeanErrorSmoothed_' phase]);
%             end
%             if saveResults
%                 entryErrorBarFig = gcf;
%                 exportfig(entryErrorBarFig,[figurename '.eps'],exportOptions)
%                 system(['epstopdf ' figurename '.eps']);
%                 system(['rm ' figurename '.eps']);
%             end
%             
%             % errorbar exit plot
%             figure;
%             shadedErrorBar(timeSeries,nanmean(allSmoothExitSpeeds,1),nanstd(allSmoothExitSpeeds,1),'r');
%             title(['mean cluster exit speeds'])
%             xlabel('frames')
%             ylabel('speed(microns/s)')
%             ylim([-150 250])
%             if enforceInClusterAfterEntryBeforeExit
%                 figurename = (['figures/entryExitSpeeds/exitSpeedsMeanErrorSmoothed_halfFiltered_' phase]);
%             else
%                 figurename = (['figures/entryExitSpeeds/exitSpeedsMeanErrorSmoothed_' phase]);
%             end
%             if saveResults
%                 exitErrorBarFig = gcf;
%                 exportfig(exitErrorBarFig,[figurename '.eps'],exportOptions)
%                 system(['epstopdf ' figurename '.eps']);
%                 system(['rm ' figurename '.eps']);
%             end
        end
    end
end