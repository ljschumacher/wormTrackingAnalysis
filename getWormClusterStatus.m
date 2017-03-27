function [ inCluster, loneWorms, rest ] = getWormClusterStatus(trajData, frame,...
    pixelsize, maxNeighbourDist, inClusterRadius, inClusterNeighbourNum)
% for a given frame, computes which worms are in/out of cluster based on
% positions
% returns logical vectors to index worms in/out of cluster

[x, y] = getWormPositions(trajData, frame);

if numel(x)>1 % need at least two worms in frame to calculate distances
    D = squareform(pdist([x y]).*pixelsize); % distance of every worm to every other
    % find lone worms
    mindist = min(D + max(max(D))*eye(size(D)));
    loneWorms = mindist>=maxNeighbourDist;
    % find worms in clusters
    numCloseNeighbours = sum(D<inClusterRadius,2);
    inCluster = numCloseNeighbours>=inClusterNeighbourNum;
    % rest is worms neither in cluster nor lone
    rest = ~loneWorms&~inCluster;
else
    inCluster = false(size(x'));
    loneWorms = false(size(x'));
    rest = true(size(x'));
end
end

