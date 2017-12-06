function continuousFrameRuns = getContinuousTraj(wormFrames)

% function takes a list of discontinuous frames belonging to a single worm
% and returns a cell array containing frames that make up continuous trajectories

%% INPUT
% wormpathFrames: vector of frame numbers (double), from a single worm_index that may have become discontinuous due to various filtering steps
%% OUTPUT
% continuousFrameRuns: cell array containing continuous trajectories

%% function

% take difference between adjacent frame numbers (continuous frames will
% give a difference of 1)
frameNumberDiff = diff(wormFrames);

% initialise
continuousFrameRuns = {}; % holds consecutive trajectories
continuousRunCtr = 1; % indicates the number of the current trajectory
continuousFrameRuns{continuousRunCtr} = wormFrames(1); % fill in the first trajectory with the first frame number value

% loop through each frame number difference value
for i = 1:numel(frameNumberDiff) 
    % if this frame continues from the previous one
    if frameNumberDiff(i) == 1 
        % add this frame to the current trajectory
        continuousFrameRuns{continuousRunCtr} = [continuousFrameRuns{continuousRunCtr} wormFrames(i+1)]; 
    else
        % start a new trajectory
        continuousRunCtr = continuousRunCtr + 1;
        % add the frame to the new trajectory
        continuousFrameRuns{continuousRunCtr} = wormFrames(i+1);
    end
end