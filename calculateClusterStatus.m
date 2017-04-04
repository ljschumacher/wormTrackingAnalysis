function [ min_neighbor_dist_rr, min_neighbor_dist_rg, num_close_neighbours_rg ] = calculateClusterStatus(trajData,trajData_g,pixelsize,inClusterRadius)
% loops through all frames and calculates number of close neighbours, and
% minimum distance to a neighbour for determining in cluster / lone status
numFrames = numel(unique(trajData.frame_number));
frames = unique(trajData.frame_number);
min_neighbor_dist_rr = cell(numFrames,1);
min_neighbor_dist_rg = cell(numFrames,1);
num_close_neighbours_rg = cell(numFrames,1);
for frameCtr = 1:numFrames
    frame = frames(frameCtr);
    %%% need to update this for 2 channel case
    [~, min_neighbor_dist_rr{frameCtr}] = ...
        getWormClusterStatus(trajData, frame, pixelsize, inClusterRadius);
    [num_close_neighbours_rg{frameCtr}, min_neighbor_dist_rg{frameCtr}] = ...
        getWormClusterStatus2channel(trajData, trajData_g, frame, pixelsize, inClusterRadius);
end
% pool data from all frames
min_neighbor_dist_rr = horzcat(min_neighbor_dist_rr{:});
min_neighbor_dist_rg = horzcat(min_neighbor_dist_rg{:});
num_close_neighbours_rg = horzcat(num_close_neighbours_rg{:});
end

