function [ neighbr_distances, min_neighbr_dist, num_close_neighbrs ] = ...
    calculateClusterStatus(trajData_a,trajData_b,pixelsize,inClusterRadius)
% loops through all frames and calculates number of close neighbrs, and
% minimum distance to a neighbr for determining in cluster / lone status
numFrames = numel(unique(trajData_a.frame_number));
frames = unique(trajData_a.frame_number);
numNeighbrDists = 10;
neighbr_distances = NaN(length(trajData_a.frame_number),numNeighbrDists);
min_neighbr_dist = NaN(size(trajData_a.frame_number));
num_close_neighbrs = zeros(size(trajData_a.frame_number));
for frameCtr = 1:numFrames
    frame = frames(frameCtr);
    frameLogIdcs = trajData_a.frame_number==frame;    % put calculated distances and neighbour numbers into the right entries
    [neighbr_distances_rr, ~, ~] = ... % use 2-channel function here to filter the red objects to which distance is calculated, but not ones distance is calculated from
        getWormClusterStatus2channel(trajData_a, trajData_a, frame, pixelsize, inClusterRadius);
    [neighbr_distances_rg, ~, ~] = ...
        getWormClusterStatus2channel(trajData_a, trajData_b, frame, pixelsize, inClusterRadius);
    % combine distances to red and green neighbours
    neighbr_distances_combined = sort([neighbr_distances_rr, neighbr_distances_rg],2);
    if size(neighbr_distances_combined,2)<numNeighbrDists
        1;
    end
    neighbr_distances(frameLogIdcs,:) = neighbr_distances_combined(:,1:numNeighbrDists);
    min_neighbr_dist(frameLogIdcs) = neighbr_distances_combined(:,1);
    num_close_neighbrs(frameLogIdcs) = sum(neighbr_distances_combined<inClusterRadius,2);
end

end

