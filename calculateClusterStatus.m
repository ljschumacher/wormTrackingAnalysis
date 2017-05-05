function [ min_neighbr_dist_rr, min_neighbr_dist_rg, num_close_neighbrs_rg ] = ...
    calculateClusterStatus(trajData,trajData_g,pixelsize,inClusterRadius)
% loops through all frames and calculates number of close neighbrs, and
% minimum distance to a neighbr for determining in cluster / lone status
numFrames = numel(unique(trajData.frame_number));
frames = unique(trajData.frame_number);
min_neighbr_dist_rr = single(size(trajData.frame_number));
min_neighbr_dist_rg = single(size(trajData.frame_number));
num_close_neighbrs_rg = uint16(size(trajData.frame_number));
for frameCtr = 1:numFrames
    frame = frames(frameCtr);
    frameLogIdcs = trajData.frame_number==frame;    % put calculated distances and neighbour numbers into the right entries
    [~, min_neighbr_dist_rr(frameLogIdcs)] = ... % use 2-channel function here to filter the red objects to which distance is calculated, but not ones distance is calculated from
        getWormClusterStatus2channel(trajData, trajData, frame, pixelsize, inClusterRadius);
    [num_close_neighbrs_rg(frameLogIdcs), min_neighbr_dist_rg(frameLogIdcs)] = ...
        getWormClusterStatus2channel(trajData, trajData_g, frame, pixelsize, inClusterRadius);
end

end

