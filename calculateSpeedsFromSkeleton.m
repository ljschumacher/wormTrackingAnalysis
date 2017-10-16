function [ speed, velocity_x, velocity_y, speedSigned ] = calculateSpeedsFromSkeleton(trajData,skelData,skelIndcs,...
    pixelsize,frameRate,filter,smoothingWindowSize)
%% calculate midbody speeds from skeletons
% INPUTS
% trajData: output from h5read(tracked_file,'/trajectorie_Data')
% skelData: output from h5read(tracked_file,'/skeleton')
% skelIndcs: which indices of the skeleton to use for calculation, eg 19:33
%        for midbody, or 1:5 for head
% pixelsize: in microns
% frameRate: in frames per second
% filter: true or false, whether to filter for discontinuous worms
% smoothing: smoothing window size (set to 0 to disable smoothing)

% centroids of midbody skeleton
x = mean(squeeze(skelData(1,skelIndcs,:)))*pixelsize;
y = mean(squeeze(skelData(2,skelIndcs,:)))*pixelsize;
% change in centroid position over time
dxdt = gradient(x)*frameRate;
dydt = gradient(y)*frameRate;
% speed and velocity
dFramedt = gradient(double(trajData.frame_number))';
speed = sqrt(dxdt.^2 + dydt.^2)./dFramedt;
velocity_x = dxdt./dFramedt;
velocity_y = dydt./dFramedt;
% signed speed calculation
% direction of segments pointing along midbody
[~, dyds] = gradient(squeeze(skelData(2,skelIndcs,:)),-1);
[~, dxds] = gradient(squeeze(skelData(1,skelIndcs,:)),-1);
% sign speed based on relative orientation of velocity to body
speedSigned = getSignedSpeed(velocity,[mean(dxds); mean(dyds)]);
if filter
% ignore first and last frames of each worm's track
wormChangeIndcs = gradient(double(trajData.worm_index_joined))~=0;
speedSigned(wormChangeIndcs)=NaN;
% ignore frames with bad skeletonization
speedSigned(trajData.is_good_skel~=1)=NaN;
% ignore skeletons otherwise filtered out
speedSigned(~trajData.filtered) = NaN;
end
if smoothingWindowSize>0
% smooth speed to denoise
speedSigned = smooth(speedSigned,smoothingWindowSize,'moving');
end

end

