function [newxcoords, newycoords] = rotateTraj(oldxcoords, oldycoords,numFramesForRotationAngle)

% function takes xy coordinates of a given trajectory as separate column
% vectors, and rotate coordinates to by the first trajectory angle

%% INPUTS:
% oldxcoords: column vector of original x-coordinates.
% oldycoords: column vector of original y-coordinates. 
% numFramesForRotationAngle: [1x1] double. The number of frames to use for
% rotation angle calculation. If unspecified then default is 3 (i.e. first three points or two segments)

%% OUTPUTS:
% newxcoords: column vector the same size as oldxoords, containing new
% x-coordinates rotated by first trajectory angle.
% newycoords: column vector the same size as oldxoords, containing new
% y-coordinates rotated by first trajectory angle.

%% FUNCTION

if nargin <3
    numFramesForRotationAngle = 3;
end

% calculate the angle given by the first trajectory segment
xdiff = oldxcoords(numFramesForRotationAngle) - oldxcoords(1);
ydiff = oldycoords(numFramesForRotationAngle) - oldycoords(1);
rotateAngle = atan2(ydiff, xdiff);

% make rotation matrix
rotationMatrix = [cos(rotateAngle), sin(rotateAngle); -sin(rotateAngle), cos(rotateAngle)];

% rotate coordinates
oldxycoords = [oldxcoords'; oldycoords']; % turn old xy coordinates into 2-row matrix for matrix multiplication
newxycoords = rotationMatrix * oldxycoords; % apply 2x2 rotation matrix to old xy coordinates
newxcoords = newxycoords(1,:);
newycoords = newxycoords(2,:);