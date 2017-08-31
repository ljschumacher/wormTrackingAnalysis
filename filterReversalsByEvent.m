function [ timeToRev, timeToRev_censored ] = ...
                filterReversalsByEvent(revStartInd, eventLogInd, worm_index, maxPostEventTime)
% filters reversals occurring after specified events, and gives back time
% to reversal after the event
% eventLogInd specifies whether event is ongoing

%find start of events
eventStartInd = find([eventLogInd(1); ~eventLogInd(1:end-1)&eventLogInd(2:end)]);
nEvents = numel(eventStartInd);
timeToRev = NaN(nEvents,1);
timeToRev_censored = false(nEvents,1);
% find the next reversal occuring after each event
for eventCtr = 1:nEvents
    thisEventInd = eventStartInd(eventCtr);
    nextRevInd = revStartInd(find(revStartInd>=thisEventInd,1,'first'));
    % check if this reversal is still from the same worm
    try
    if worm_index(thisEventInd)==worm_index(nextRevInd)
        timeToRev(eventCtr) = nextRevInd - thisEventInd;
    else % this means the track identity was lost before a reversal occurred
        if isempty(nextRevInd), nextRevInd = numel(worm_index); end
        % find the last frame of the track of the same worm after the event
        lastTracked = find(worm_index(1:nextRevInd)==worm_index(thisEventInd),1,'last');
        % keep time until track was lost, but mark it as right-censored
        timeToRev(eventCtr) = lastTracked - thisEventInd;
        timeToRev_censored(eventCtr) = true;
    end
    if timeToRev(eventCtr)>maxPostEventTime
        timeToRev(eventCtr) = maxPostEventTime;
        timeToRev_censored(eventCtr) = true;
    end
    catch
        1;
    end
end
assert(~any(timeToRev<0),'Error: negative times until reversals after event')