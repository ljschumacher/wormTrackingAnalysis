function [ min_neighbr_dist_rr, min_neighbr_dist_rg, num_close_neighbrs_rg ] = ...
    calculateClusterStatus(trajData,trajData_g,pixelsize,inClusterRadius)
% loops through all frames and calculates number of close neighbrs, and
% minimum distance to a neighbr for determining in cluster / lone status
numFrames = numel(unique(trajData.frame_number));
frames = unique(trajData.frame_number);
min_neighbr_dist_rr = cell(numFrames,1);
min_neighbr_dist_rg = cell(numFrames,1);
num_close_neighbrs_rg = cell(numFrames,1);
for frameCtr = 1:numFrames
    frame = frames(frameCtr);
    [~, min_neighbr_dist_rr{frameCtr}] = ... % use 2-channel function here to filter the red objects to which distance is calculated, but not ones distance is calculated from
        getWormClusterStatus2channel(trajData, trajData, frame, pixelsize, inClusterRadius);
    [num_close_neighbrs_rg{frameCtr}, min_neighbr_dist_rg{frameCtr}] = ...
        getWormClusterStatus2channel(trajData, trajData_g, frame, pixelsize, inClusterRadius);
end
% pool data from all frames
min_neighbr_dist_rr = horzcat(min_neighbr_dist_rr{:});
min_neighbr_dist_rg = horzcat(min_neighbr_dist_rg{:});
num_close_neighbrs_rg = horzcat(num_close_neighbrs_rg{:});
end

