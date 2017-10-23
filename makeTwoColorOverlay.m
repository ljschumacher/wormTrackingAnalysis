% script takes hdf5 masked movies, reads them frame by frame, and overlay the two
% channels to generate new frame and avi movie.

recordingNumber = 51.2;
startFrame = 12500;
endFrame = 24250;
frameInterval = 1;

filename_g = rdir(['*/recording' num2str(recordingNumber) 'g*/*.hdf5']);
filename_r = rdir(['*/recording' num2str(recordingNumber) 'r*/*.hdf5']);

newfilename = ['recording' recordingNumber 'rg_f' num2str(startFrame) '-f' num2str(endFrame) '.avi'];

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
open(vid_rg)

% loop through each frame
for frameNumber = startFrame:frameInterval:endFrame

    % load the frame (frames are greyscale, equally distributed in all RGB channels)
    frame_g = h5read(filename_g, '/mask', [1, 1, frameNumber], [dims_g(1), dims_g(2), 1]);
    frame_r = h5read(filename_r, '/mask', [1, 1, frameNumber], [dims_r(1), dims_r(2), 1]);
    % apply offset
    frame_g = frame_g(yOffset+1:end,1:end-xOffset,:);
    frame_r = frame_r(1:end-yOffset,xOffset+1:end,:);
    % add frame number/time stamp to green image
    frame_g = AddTextToImage(frame_g,['frame ' num2str(frameNumber)],[100 100],[255,255,255],'Arial',50);
    % make new overlay frame
    frame_rg = cat(3,frame_r(:,:,1),frame_g(:,:,1),frame_g(:,:,1));
    % write new frame to video
    writeVideo(vid_rg, frame_rg)
end

close(vid_rg)