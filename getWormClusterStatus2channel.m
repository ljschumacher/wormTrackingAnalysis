function [ neighbr_distances, numCloseNeighbrs, mindist ] = getWormClusterStatus2channel(...
    trajData_a,trajData_b, frame,...
    pixelsize, inClusterRadius)
% for a given frame, computes which worms are in/out of cluster based on
% positions
% returns logical vectors to index worms in/out of cluster
% 2 channel version, calculates cluster status for red worms based on
% distances to green worms
% currently ignores neighbourship of red worms to each other
% issues / todo:
% - could also use knnsearch for nearest neighbours

[x_a, y_a] = getWormPositions(trajData_a, frame, false);% don't use filtered objects, as we want as many data as entries in the original
[x_b, y_b] = getWormPositions(trajData_b, frame, true); % do filter for the object to which the distance is calculated (don't care if close to some artefact)

if numel(x_b)>=1&&numel(x_a)>=1 % need at least two worms in frame to calculate distances
    a2bDistances = pdist2([x_a y_a],[x_b y_b]).*pixelsize; % distance of every red worm to every green
    % sometimes we may be using this to calculate distances where a and b
    % are the same. for this case set distances to self to Inf
    a2bDistances(a2bDistances==0) = Inf;
    % find lone worms
    neighbr_distances = sort(a2bDistances,2);
    numNeighbrs = size(neighbr_distances,2);
    if numNeighbrs>=10
        neighbr_distances = neighbr_distances(:,1:10);
    else
        neighbr_distances(:,numNeighbrs+1:10) = NaN;
    end
    mindist = neighbr_distances(:,1)';  % transpose for consistency with getWormClusterStatus.m
    % find worms in clusters
    numCloseNeighbrs = sum(a2bDistances<inClusterRadius,2)'; % transpose for conistency with getWormClusterStatus.m
else
    numCloseNeighbrs = zeros(size(x_a'));
    mindist = NaN(size(x_a'));
    neighbr_distances = NaN(length(x_a),10);
end
end

