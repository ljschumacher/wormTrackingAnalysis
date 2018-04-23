function [ x, y] = getSkelPositions(skelData, currentFrameLogInd)
% returns the x, y positions of skeleta for a given frame

% get positions
x = squeeze(skelData(1,:,currentFrameLogInd));
y = squeeze(skelData(2,:,currentFrameLogInd));
% reshape to have expected dimensions N by M
M = size(skelData,2);
N = nnz(currentFrameLogInd);
if size(x,1)==M&&size(x,2)==N
    x = x';
    y = y';
end
end

