function [ numCloseNeighbrs, mindist] = getWormClusterStatus(trajData, frame,...
    pixelsize, inClusterRadius)
% for a given frame, computes which worms are in/out of cluster based on
% positions
% returns logical vectors to index worms in/out of cluster

[x, y] = getWormPositions(trajData, frame, false);

if numel(x)>1 % need at least two worms in frame to calculate distances
    D = squareform(pdist([double(x) double(y)]).*double(pixelsize)); % distance of every worm to every other
    % find lone worms
    mindist = min(D + max(max(D))*eye(size(D)));
    % find worms in clusters
    numCloseNeighbrs = sum(D<inClusterRadius,2);
else
    numCloseNeighbrs = zeros(size(x'));
    mindist = NaN(size(x'));
end
end

