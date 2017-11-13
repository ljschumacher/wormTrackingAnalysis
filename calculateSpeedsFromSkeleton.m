function [ speed, velocity_x, velocity_y, speedSigned, acceleration_x, acceleration_y ] =...
    calculateSpeedsFromSkeleton(trajData,skelData,skelIndcs,...
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
dxdFrame = gradient(x)*frameRate;
dydFrame = gradient(y)*frameRate;
% speed and velocity
dFramedt = gradient(double(trajData.frame_number))';
speed = sqrt(dxdFrame.^2 + dydFrame.^2)./dFramedt;
velocity_x = dxdFrame./dFramedt;
velocity_y = dydFrame./dFramedt;
% acceleration - change in velocity
acceleration_x = gradient(velocity_x)./dFramedt;
acceleration_y = gradient(velocity_y)./dFramedt;
% signed speed calculation
% direction of segments pointing along midbody
[~, dyds] = gradient(squeeze(skelData(2,skelIndcs,:)),-1);
[~, dxds] = gradient(squeeze(skelData(1,skelIndcs,:)),-1);
% sign speed based on relative orientation of velocity to body
speedSigned = getSignedSpeed([velocity_x; velocity_y],[mean(dxds); mean(dyds)]);
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
    speedSigned = smoothdata(speedSigned,'movmean',smoothingWindowSize,'omitnan');
end

end

