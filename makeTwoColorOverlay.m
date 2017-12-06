% script takes hdf5 masked movies, reads them frame by frame, and overlay the two
% channels to generate new frame and avi movie.

addpath('auxiliary/')
%% set parameters
frameInterval = 1;

%% load manual annotation file and extract recording and frame annotations
[~,~,annotations] = xlsread(['datalists/entryExitEvents_npr1_40_joining.xlsx'],1,'A2:I83','basic');
recordingNumbers = annotations(:,1);

%% loop through each entry and exit event
for eventCtr = 1:length(recordingNumbers)
    
    % get movie number
    recordingNumber = recordingNumbers{eventCtr};
    recordingNumber = recordingNumber(2:end); % to remove the 'r' character at the start
    
    % get start and end frame numbers
    startFrame = annotations{eventCtr,7}-200; % go 200 frames before hand annotated event start
    endFrame = annotations{eventCtr,8}+200; % go 200 frames after hand annotated event finish
    
    
    % search for the masked video file
    filelist_g = rdir(['/data2/shared/data/twoColour/MaskedVideos/*/recording' recordingNumber 'g*/*.hdf5']);
    filelist_r = rdir(['/data2/shared/data/twoColour/MaskedVideos/*/recording' recordingNumber 'r*/*.hdf5']);
    filename_g = filelist_g.name;
    filename_r = filelist_r.name;
    
    % name the new file
    newfilename = ['/data2/shared/data/twoColour/entryExitExerptsWithSpeedGraph/recording' recordingNumber 'rg_f' num2str(startFrame) '-f' num2str(endFrame) '.avi'];
    
    % allow offset to align misaligned channels
    xOffset = 25;
    yOffset = 25;
    
    % get the dimensions of the video
    fileInfo_g = hdf5info(filename_g);
    fileInfo_r = hdf5info(filename_r);
    dims_g = fileInfo_g.GroupHierarchy.Datasets(2).Dims;
    dims_r = fileInfo_r.GroupHierarchy.Datasets(2).Dims;
    
    % create output video
    vid_rg = VideoWriter(newfilename,'Motion JPEG AVI');
    vid_rg.FrameRate = 27; % speed up 3x compared to 9fps acquisition 
    open(vid_rg)
    
    % loop through each frame
    for frameNumber = startFrame:frameInterval:endFrame
        
        % load the frame (frames are greyscale, equally distributed in all RGB channels)
        frame_g = h5read(filename_g, '/mask', [1, 1, frameNumber], [dims_g(1), dims_g(2), 1]);
        frame_r = h5read(filename_r, '/mask', [1, 1, frameNumber], [dims_r(1), dims_r(2), 1]);
        % apply offset and rotate
        frame_g = frame_g(yOffset+1:end,1:end-xOffset,:)';
        frame_r = frame_r(1:end-yOffset,xOffset+1:end,:)';
        % add frame number/time stamp to green image
        frame_g = AddTextToImage(frame_g,['frame ' num2str(frameNumber)],[100 100],[255,255,255],'Arial',50);
        % make new overlay frame
        frame_rg = cat(3,frame_r(:,:,1),frame_g(:,:,1),frame_g(:,:,1));

        % write new frame to video
        writeVideo(vid_rg, frame_rg)
    end
    
    close(vid_rg)
    
end
