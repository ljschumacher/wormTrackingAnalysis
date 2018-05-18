function [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr)

% function takes pre-specified phaseFrame matrix, phase, and file counter,
% and returns first and last frame numbers for the phase of interest. 

%% INPUTS: 
% phaseFrames: matrix of phaseFrames as per manual annotation in the .xlsx datalist file. First and second columns denote the start and 
% end of the joining phase, the third and fourth columns denote the start and end of the sweeping phase. Each row corresponds to a different
% file replicate. 
% phase: string to specify 'fullMovie','joining' or 'sweeping'.
% fileCtr: [1x1] double, specifies the position of the current file of interest on the list;. Usually 1-12

%% OUTPUTS:
% firstFrame: [1x1] double,frame number for the beginning of the specified phase.
% lastFrame: [1x1] double, frame number for the end of the specified phase.

%% FUNCTION:
% correct for python indexing at 0 (i.e. in the tracking data .hdf5 file each movie starts at frame 0, rather than 1)
phaseFrames = phaseFrames-1; 

% obtain desired first and last frames accordingly
if strcmp(phase, 'fullMovie')
    firstFrame = 0;
    lastFrame = phaseFrames(fileCtr,4);
elseif strcmp(phase,'joining')
    firstFrame = phaseFrames(fileCtr,1);
    lastFrame = phaseFrames(fileCtr,2);
elseif strcmp(phase,'sweeping')
    firstFrame = phaseFrames(fileCtr,3);
    lastFrame = phaseFrames(fileCtr,4);
elseif strcmp(phase,'nontransient')
    firstFrame = 0;
    lastFrame = phaseFrames(fileCtr,4);
end