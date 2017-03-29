function [ numCloseNeighbours, mindist ] = getWormClusterStatus2channel(...
    trajData_r,trajData_g, frame,...
    pixelsize, inClusterRadius)
% for a given frame, computes which worms are in/out of cluster based on
% positions
% returns logical vectors to index worms in/out of cluster
% 2 channel version, calculates cluster status for red worms based on
% distances to green worms
% currently ignores neighbourship of red worms to each other

[x_r, y_r] = getWormPositions(trajData_r, frame);
[x_g, y_g] = getWormPositions(trajData_g, frame);

if numel(x_g)>=1&&numel(x_r)>=1 % need at least two worms in frame to calculate distances
    redToGreenDistances = pdist2([x_r y_r],[x_g y_g]).*pixelsize; % distance of every red worm to every green
    % find lone worms
    mindist = min(redToGreenDistances,[],2)'; % transpose for conistency with getWormClusterStatus.m
    % find worms in clusters
    numCloseNeighbours = sum(redToGreenDistances<inClusterRadius,2)'; % transpose for conistency with getWormClusterStatus.m
else
    numCloseNeighbours = zeros(size(x_r'));
    mindist = NaN(size(x_r'));
end
end

