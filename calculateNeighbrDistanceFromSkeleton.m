function [ neighbr_distances_skeleton ] = ...
    calculateNeighbrDistanceFromSkeleton(skelData_r,trajData_r,trajData_g,pixelsize)
% loops through all frames and calculates distance to nearest ten
% neighbours for each part of the skeleton
frames = unique(trajData_r.frame_number);
numFrames = numel(frames);
numNeighbrDists = 10;
M = size(skelData_r,2);
neighbr_distances_skeleton = NaN(numFrames,M,numNeighbrDists);
for frameCtr = 1:numFrames
    frame = frames(frameCtr);
    frameLogIdcs = trajData_r.frame_number==frame;    % put calculated distances into the right entries
    neighbr_distances_rr = getWormNbrDistancesSkeleton(skelData_r, trajData_r, frame,frameLogIdcs, pixelsize);
    if ~isempty(trajData_g)
        neighbr_distances_rg = getWormNbrDistancesSkeleton(skelData_r, trajData_g, frame, frameLogIdcs, pixelsize);
        % combine distances to red and green neighbours
        neighbr_distances_combined = sort(cat(3,neighbr_distances_rr, neighbr_distances_rg),3);
    else
        neighbr_distances_combined = neighbr_distances_rr;
    end
    neighbr_distances_skeleton(frameLogIdcs,:,:) = neighbr_distances_combined(:,:,1:numNeighbrDists);
end

end

