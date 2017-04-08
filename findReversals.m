function [revStartInd, revDuration] = findReversals(midbody_speed, worm_index)
% finds reversals based on worm tracking data features file

% define reversal starting events as when sign of speed changes from + to -
revStartInd = find(midbody_speed(1:end-1)>0&midbody_speed(2:end)<0);
% reversals end when sign of speed changes from - to + (but there could
% be more or fewer of these than reversal start times)
revEndInd = 1+find(midbody_speed(1:end-1)<0&midbody_speed(2:end)>0);
% match up reversal end times with starts
revEndIndMatched = NaN(size(revStartInd));
nRevs = length(revStartInd);
for revCtr = 1:nRevs
    if revStartInd(revCtr)<max(revEndInd)
        revEndIndMatched(revCtr) = revEndInd(find(revEndInd>revStartInd(revCtr),1,'first'));
        if worm_index(revStartInd(revCtr)) ~= worm_index(revEndIndMatched(revCtr))
            % reversal 'end' is actually from a different worm, which means
            % that the reversal end is not tracked (at least not under the
            % same identity)
            revEndIndMatched(revCtr) = NaN;
        end
        % it can happen that the end of the reversal is not tracked, and
        % that another reversal has started before the matched end
        if revCtr<nRevs&&revEndIndMatched(revCtr)>revStartInd(revCtr+1)
            revEndIndMatched(revCtr) = NaN;
        end 
    else % reversal still ongoing when tracking ends
        revEndIndMatched(revCtr) = NaN;
    end
end
revDuration = revEndIndMatched - revStartInd;
end

