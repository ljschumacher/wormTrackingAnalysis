function [ x, y] = getWormPositions(skelTrajData, framenumber, filter)
% calculates worm speeds from tracking data for a given frame
% returns the x, y positions and u, v displacement components
if nargin<3
    filter = false;
end
currentFrameLogInd = skelTrajData.frame_number==framenumber;
% optionally filter blobs
if isfield(skelTrajData,'filtered')&&filter
    currentFrameLogInd = currentFrameLogInd&skelTrajData.filtered;
end
% get positions
x = skelTrajData.coord_x(currentFrameLogInd);
y = skelTrajData.coord_y(currentFrameLogInd);
end

