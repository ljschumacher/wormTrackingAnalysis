function [ timeToRev, timeToRev_censored, preRevEventInd ] = ...
    filterReversalsByEvent(revStartInd, eventLogInd, worm_index, maxPostEventTime)
% filters reversals occurring after specified events, and gives back time
% to reversal after the event
% eventLogInd specifies whether event is ongoing

% issues:
% this will not count the times when an event did not result in a reversal.
% to do that, one would need to loop through events, but be careful never
% to count a reversal twice for occuring after different events

%find start of events
eventStartInd = find([eventLogInd(1); ~eventLogInd(1:end-1)&eventLogInd(2:end)]);
% loop through reversals, and find the last preceding event, making sure it's from the same worm
nReversals = numel(revStartInd);
% pre-allocate more than enough spave
timeToRev = NaN(nReversals,1);
timeToRev_censored = false(nReversals,1);
preRevEventInd = NaN(nReversals,1);
for revCtr = 1:nReversals
    thisRevInd = revStartInd(revCtr);
    prevEventInd = eventStartInd(find(eventStartInd< thisRevInd,1,'last'));
    % check event and reversal are from the same worm
    if worm_index(thisRevInd) == worm_index(prevEventInd)
        timeToRev(revCtr) = thisRevInd - prevEventInd;
        preRevEventInd(revCtr) = prevEventInd;
        % if the track identity was lost before a reversal occurred, we disregard it
        if timeToRev(revCtr)>maxPostEventTime % count as track did not end in a reversal
            timeToRev(revCtr) = maxPostEventTime;
            timeToRev_censored(revCtr) = true;
        end
    end
end
% now keep only those reversals with preceding events
revKeepLogInd = ~isnan(timeToRev);
timeToRev = timeToRev(revKeepLogInd);
timeToRev_censored = timeToRev_censored(revKeepLogInd);
preRevEventInd = preRevEventInd(revKeepLogInd);

%     %% old code % THIS CAN LEAD TO REVERSALS BEING ATTRIBUTED TO MULTIPLE EVENTS!
%     nEvents = numel(eventStartInd);
% timeToRev = NaN(nEvents,1);
% timeToRev_censored = false(nEvents,1);
% for eventCtr = 1:nEvents
%     thisEventInd = eventStartInd(eventCtr);
%     nextRevInd = revStartInd(find(revStartInd>=thisEventInd,1,'first'));
%     % check if this reversal is still from the same worm
%     if worm_index(thisEventInd)==worm_index(nextRevInd)
%         timeToRev(eventCtr) = nextRevInd - thisEventInd;
%     else % this means the track identity was lost before a reversal occurred
%         if isempty(nextRevInd), nextRevInd = numel(worm_index); end
%         % find the last frame of the track of the same worm after the event
%         lastTracked = find(worm_index(1:nextRevInd)==worm_index(thisEventInd),1,'last');
%         % keep time until track was lost, but mark it as right-censored
%         timeToRev(eventCtr) = lastTracked - thisEventInd;
%         timeToRev_censored(eventCtr) = true;
%     end
%     if timeToRev(eventCtr)>maxPostEventTime % count as track did not end in a reversal
%         timeToRev(eventCtr) = maxPostEventTime;
%         timeToRev_censored(eventCtr) = true;
%     end
% end
assert(~any(timeToRev<0),'Error: negative times until reversals after event')