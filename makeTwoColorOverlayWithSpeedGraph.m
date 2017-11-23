% script takes hdf5 masked movies, reads them frame by frame, and overlay the two
% channels to generate new frame and avi movie.

%% set parameters
preExitDuration = 20; % duration (in seconds) before a worm exits a cluster to be included in the leave cluster analysis
postExitDuration = 20; % duration (in seconds) after a worm exits a cluster to be included in the leave cluster analysis
frameRate = 9; % frameRate is 9 for this twoColour dataset
manualEventMaxDuration = 400; % max number of frames that contains the beginning and end of an annotated event
midbodyIndcs = 19:33;
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
smoothWindow = 9; % number of frames to smooth over for later trajectory-specific midbody speed calculations

exportOptions = struct('Format','EPS2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',20,...
    'LineWidth',3);

%% load manual annotation files

% entry exit event annotations
[~,~,annotations] = xlsread(['/data2/shared/data/twoColour/entryExitEvents_npr1_40_joining.xlsx'],1,'A2:I80','basic');
% frames for movie phase of interest, for extract joining phase limits
[phaseFrames,filenames,~] = xlsread(['datalists/npr1_40_r_list.xlsx'],1,'A1:E15','basic');

% trim annotations down to those of just entry and exit events (i.e. ignore "roaming past" cases etc.)
entryExitEventLogInd = false(1,size(annotations,1));
for eventCtr = 1:size(annotations,1)
    if strcmp(annotations{eventCtr,5},'enter') | strcmp(annotations{eventCtr,5},'exit')
        entryExitEventLogInd(eventCtr) = true;
    end
end
annotations = annotations(entryExitEventLogInd,:);
recordingNumbers = annotations(:,1);

%% loop through each entry and exit event

% generate timeSeries (x) for speed plot
timeSeries = [-preExitDuration*frameRate:postExitDuration*frameRate+manualEventMaxDuration]/frameRate;

for eventCtr = 1%:length(recordingNumbers)
    
    %% get annotation info and find files
    
    % get movie number
    recordingNumber = recordingNumbers{eventCtr};
    recordingNumber = recordingNumber(2:end); % to remove the 'r' character at the start
    
    firstPhaseFrame = 0;
    lastPhaseFrame = 32399;
    
%     % get joining phase frames for this movie
%     for fileCtr = 1:length(filenames)
%         filename = filenames(fileCtr);
%         if ~isempty(strfind(filename{1},recordingNumber))
%             firstPhaseFrame = phaseFrames(fileCtr,1);
%             lastPhaseFrame = phaseFrames(fileCtr,2);
%         end
%     end
    
    % get worm index, plus start and end frame numbers
    wormIndex = annotations{eventCtr,2};
    startFrame = annotations{eventCtr,7}-preExitDuration*frameRate; % go 20s before hand annotated event start
    endFrame = annotations{eventCtr,8}+postExitDuration*frameRate; % go 20s after hand annotated event finish
    
    % search for the masked video and skeletons file
    filelist_MaskedVideo_g = rdir(['/data2/shared/data/twoColour/MaskedVideos/*/recording' recordingNumber 'g*/*.hdf5']);
    filelist_MaskedVideo_r = rdir(['/data2/shared/data/twoColour/MaskedVideos/*/recording' recordingNumber 'r*/*.hdf5']);
    filename_MaskedVideo_g = filelist_MaskedVideo_g.name;
    filename_MaskedVideo_r = filelist_MaskedVideo_r.name;
    filelist_Skeletons_r = rdir(['/data2/shared/data/twoColour/Results/*/recording' recordingNumber 'r*/*_X1_skeletons.hdf5']);
    filename_Skeletons_r = filelist_Skeletons_r.name;
    
    %% calculate midbody speed from manually joined red skeletons file
    
    % load files
    trajData = h5read(filename_Skeletons_r,'/trajectories_data');
    skelData = h5read(filename_Skeletons_r,'/skeleton'); % in pixels

    % midbody speed (y) calculation
    % load unsorted xy coords; sort later using worm index for interpolation
    xcoords = squeeze(skelData(1,:,:));
    ycoords = squeeze(skelData(2,:,:));
    
    % skip data filtering
    
    % trim start and end frames
    % check that the extended start doesn't go beyond the phase of interest
    if startFrame <= firstPhaseFrame
        % take note of omitted frames for alignment purposes
        beforeStartFrameNum = firstPhaseFrame-startFrame;
        % exclude frames before the start of the specified phase
        startFrame = firstPhaseFrame;
    end
    % check that the extended start doesn't go beyond the phase of interest
    if  endFrame > lastPhaseFrame
        % take note of omitted frames for alignment purposes
        afterEndFrameNum =  endFrame -lastPhaseFrame;
        % exclude frames beyond the end of the specified phase
        endFrame = lastPhaseFrame;
    end
    
    % get aligned list of frames for the event
    thisEventSpeeds = NaN(1,frameRate*(preExitDuration+postExitDuration)+1+manualEventMaxDuration);
    thisEventFrames = NaN(1,frameRate*(preExitDuration+postExitDuration)+1+manualEventMaxDuration);
    if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
        thisEventFrames(beforeStartFrameNum+1:end) = startFrame:endFrame;
    elseif exist('afterEndFrameNum','var')
        thisEventFrames(1:end-afterEndFrameNum) = startFrame:endFrame;
    else
        thisEventFrames = startFrame:endFrame;
    end
    
    % sort midbody speed and frames for the worm at consideration
    % get the indices for the worm of interest
    uniqueWormIdx = find(trajData.worm_index_manual == wormIndex);
    % extract skeleton data for worm of interest
    xcoordsSorted = xcoords(:,uniqueWormIdx);
    ycoordsSorted = ycoords(:,uniqueWormIdx);
    frameNumberSorted = trajData.frame_number(uniqueWormIdx);
    
    % interpolate over NaN values for sorted xy coordinates
    for nodeCtr = 1:size(xcoords,1)
        xcoordsNode = xcoordsSorted(nodeCtr,:);
        ycoordsNode = ycoordsSorted(nodeCtr,:);
        xcoordsNode = naninterp(xcoordsNode); % naninterp only works for vectors so go node by node
        ycoordsNode = naninterp(ycoordsNode);
        xcoordsSorted(nodeCtr,:) = xcoordsNode;
        ycoordsSorted(nodeCtr,:) = ycoordsNode;
    end
    
    % calculate midbodyspeed using sorted, interpolated xy coordinates
    % centroids of midbody skeleton
    x = mean(xcoordsSorted(midbodyIndcs,:))*pixelsize;
    y = mean(ycoordsSorted(midbodyIndcs,:))*pixelsize;
    % change in centroid position over time
    dxdt = gradient(x)*frameRate;
    dydt = gradient(y)*frameRate;
    % speed
    dFramedt = gradient(double(frameNumberSorted))';
    midbodySpeed = sqrt(dxdt.^2 + dydt.^2)./dFramedt;
    
    % go through each frame to obtain speed
    for frameCtr = 1:length(thisEventFrames)
        frameNumber = thisEventFrames(frameCtr);
        if ~isnan(frameNumber)
            wormFrameLogInd = frameNumberSorted == frameNumber;
            if nnz(wormFrameLogInd)~=0
                assert(nnz(wormFrameLogInd) ==1);
                thisEventSpeeds(frameCtr) = midbodySpeed(wormFrameLogInd);% obtain speed
            end
        end
    end
    
    % set maximum speed
    thisEventSpeeds(thisEventSpeeds>1500) = NaN;
    
    % smooth speeds
    thisEventSpeeds = smoothdata(thisEventSpeeds,2,'movmean',smoothWindow); %%%% check syntax
    
    % clear variables
    clear beforeStartFrameNum
    clear afterEndFrameNum
    
    %% create output file
    
    % name the new file
    newfilename = ['recording' recordingNumber 'rggraph_f' num2str(startFrame) '-f' num2str(endFrame) '.avi'];
    
    % allow offset to align misaligned channels
    xOffset = 25;
    yOffset = 25;
    
    % get the dimensions of the video
    fileInfo_g = hdf5info(filename_MaskedVideo_g);
    fileInfo_r = hdf5info(filename_MaskedVideo_r);
    dims_g = fileInfo_g.GroupHierarchy.Datasets(2).Dims;
    dims_r = fileInfo_r.GroupHierarchy.Datasets(2).Dims;
    
    % create output video
    vid_rg = VideoWriter(newfilename,'Motion JPEG AVI');
    vid_rg.FrameRate = 27; % speed up 3x compared to 9fps acquisition 
    open(vid_rg)
    
    % loop through each frame
    for frameCtr = 1:length(thisEventFrames)
        frameNumber = thisEventFrames(frameCtr);
        % load the frame (frames are greyscale, equally distributed in all RGB channels)
        frame_g = h5read(filename_MaskedVideo_g, '/mask', [1, 1, frameNumber], [dims_g(1), dims_g(2), 1]);
        frame_r = h5read(filename_MaskedVideo_r, '/mask', [1, 1, frameNumber], [dims_r(1), dims_r(2), 1]);
        % apply offset and rotate
        frame_g = frame_g(yOffset+1:end,1:end-xOffset,:)';
        frame_r = frame_r(1:end-yOffset,xOffset+1:end,:)';
        % add time stamp to green image
        frame_g = AddTextToImage(frame_g,['time = ' sprintf('%0.2f',frameNumber/9) 's'],[100 100],[255,255,255],'Arial',50);
        % make new overlay frame
        frame_rg = cat(3,frame_r(:,:,1),frame_g(:,:,1),frame_g(:,:,1));
        
        %%%%%%%%%% overlay frame_rg with speed graph somehow...%%%%%%%%%%%%
        
        
        % plot speed graph 
        speedFig = figure; hold on 
        plot(timeSeries,thisEventSpeeds,'r')
        plot(timeSeries(frameCtr),thisEventSpeeds(frameCtr),'bo')
        vline(0,'k')
        title(['event ' num2str(eventCtr)])
        xlabel('time(s)')
        ylabel('speed(microns/s)')
        xlim([-preExitDuration postExitDuration])
        ylim([0 450])
        
        % write new frame to video
        writeVideo(vid_rg, frame_rg)
    end
    
    close(vid_rg)
    
end
