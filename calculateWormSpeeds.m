function [ x, y, u, v ] = calculateWormSpeeds(trajData, framenumber, filter)
% calculates worm speeds from tracking data for a given frame
% returns the x, y positions and u, v displacement components
if nargin<3
    filter = false;
end
currentFrameLogInd = trajData.frame_number==framenumber;
nextFrameLogInd = trajData.frame_number==framenumber+1;
% optionally filter blobs
if filter&&isfield(trajData,'filtered')
    currentFrameLogInd = currentFrameLogInd&trajData.filtered;
    nextFrameLogInd = nextFrameLogInd&trajData.filtered;
end
wormInds = unique(trajData.worm_index_joined(currentFrameLogInd));
nWorms = numel(wormInds);
% get positions
x = trajData.coord_x(currentFrameLogInd);
y = trajData.coord_y(currentFrameLogInd);
% initialise displacements
u = NaN(size(x));
v = NaN(size(y));
% calculate displacements, checking that there are no tracking errors
if nnz(currentFrameLogInd)>nWorms
    warning(['Some worm(s) appear(s) more than once in frame ' num2str(framenumber) '. Cannot calculate speed.'])
else
    % check which entries are from the next frame and contain current worms
    nextFrameWormLogInd = nextFrameLogInd&ismember(trajData.worm_index_joined,wormInds);
    % check which of the current worms are in the next frame
    overlapLogInd = ismember(wormInds,trajData.worm_index_joined(nextFrameLogInd));
    if nnz(nextFrameWormLogInd) > nnz(overlapLogInd)
        warning(['Some worm(s) from ' num2str(framenumber) ' appear(s) more than once in frame ' num2str(framenumber) '+1. Cannot calculate speed.'])
    else
        u(overlapLogInd) = trajData.coord_x(nextFrameWormLogInd) - x(overlapLogInd);
        v(overlapLogInd) = trajData.coord_y(nextFrameWormLogInd) - y(overlapLogInd);
        if any(~overlapLogInd)
            disp([ num2str(nnz(~overlapLogInd)) ' of ' num2str(nWorms) ' worms did not appear in frame ' num2str(framenumber) '+1.'])
        end
    end
end
end

