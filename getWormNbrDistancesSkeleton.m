function [ neighbr_distances ] = getWormNbrDistancesSkeleton(...
    skelData,trajData, frame,frameLogIdcs,pixelsize)
% for a given frame, computes distances to nearest nbrs from each point in
% skeleton
M = size(skelData,2);
[x_skel, y_skel] = getSkelPositions(skelData, frameLogIdcs);% don't use filtered objects, as we want as many data as entries in the original
[x_blobs, y_blobs] = getWormPositions(trajData, frame, true); % do filter for the object to which the distance is calculated (don't care if close to some artefact)

neighbr_distances = NaN(size(x_skel,1),M,10);
    
if numel(x_blobs)>=1&&numel(x_skel)>=1 % need at least two worms in frame to calculate distances
    for nodeCtr = 1:M
        a2bDistances = pdist2([x_skel(:,nodeCtr) y_skel(:,nodeCtr)],[x_blobs y_blobs]).*pixelsize; % distance of every skeleton node to every blob
        sorted_distances = sort(a2bDistances,2);
        numNeighbrs = size(sorted_distances,2);
        if numNeighbrs>=10
            neighbr_distances(:,nodeCtr,:) = sorted_distances(:,1:10);
        else
            neighbr_distances(:,nodeCtr,1:numNeighbrs) = sorted_distances;
        end
    end
end
end

