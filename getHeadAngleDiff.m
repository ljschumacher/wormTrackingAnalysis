function [headAngleDiff,framesElapsed] = getHeadAngleDiff(xcoords, ycoords, marker, smoothing, frameRate)

% function returns a vector of head angle differences (in radian) between
% each frame of a given trajectory

%% INPUTS:
% xcoords: vector or matrix (if multiple nodes) containing x coordinates for a given trajectory.
% ycoords: vector or matrix (if multiple nodes) containing y coordinates for a given trajectory.
% marker: string, 'pharynx' or 'bodywall'.
% smoothing: logical, true or false.
% frameRate: [1x1] scalar, 3 or 9.
%% OUTPUTS:
% headAngleDiff: vector of head angle differences (in radian) between each frame of the given trajectory
% framesElapsed: [1x1] scalar indicating the number of frames of the trajectory

%% FUNCTION

if nargin<4
    smoothing = false;
    frameRate = [];
    if nargin<3
        marker = 'notBodyWall';
    end
end

% calculate head angles
[angleArray,meanAngles] = makeAngleArray(xcoords,ycoords);
headAngle = angleArray+meanAngles;
% take mean head angles for body wall marker
if strcmp(marker,'bodywall')
    headAngle = nanmean(headAngle,2); % take mean of the 8 head nodes
end
% get angle difference
if smoothing
    smoothFactor = frameRate+1; % set smoothing to be over 1 second
    angleDiff = headAngle(2:end) - headAngle(1:end-1);
    totalSmoothedFrames = length(angleDiff)-smoothFactor;
    headAngleDiff = NaN(totalSmoothedFrames,1);
    for smoothFrameCtr = 1:totalSmoothedFrames
        headAngleDiff(smoothFrameCtr) = nanmean(angleDiff(smoothFrameCtr:smoothFrameCtr+smoothFactor));
    end
    framesElapsed = length(angleDiff);
else
    headAngleDiff = headAngle(2:end) - headAngle(1:end-1);
    framesElapsed = length(headAngleDiff);
end
