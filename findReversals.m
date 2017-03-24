function [revStartInd, revDuration] = findReversals(features)
% finds reversals based on worm tracking data features file

% define reversal starting events as when sign of speed changes from + to -
revStartInd = find(features.midbody_speed(1:end-1)>0&features.midbody_speed(2:end)<0);
% reversals end when sign of speed changes from - to + -- but there could
% be more of these than reversal start times
revEndInd = 1+find(features.midbody_speed(1:end-1)<0&features.midbody_speed(2:end)>0);
% match up reversal end times with starts
revEndIndMatched = NaN(size(revStartInd));
for revCtr = 1:length(revStartInd)
    if revStartInd(revCtr)<max(revEndInd)
        revEndIndMatched(revCtr) = revEndInd(find(revEndInd>revStartInd(revCtr),1,'first'));
        if features.worm_index(revStartInd(revCtr)) ~= features.worm_index(revEndIndMatched(revCtr))
            % reversal 'end' is actually from a different worm, which means
            % that the reversal end is not tracked (at least not under the
            % same identity)
            revEndIndMatched(revCtr) = NaN;
        end
    else % reversal still ongoing when tracking ends
        revEndIndMatched(revCtr) = NaN;
    end
end
revDuration = revEndIndMatched - revStartInd;
end

