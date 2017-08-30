% Script loads saved xy coordinate data and plots a number of sample
% trajectories starting at the origin and rotated. 
% Script gives the option to apply a minimum trajectory length filter.

clear
% close all

% set parameters
dataset = 2;
wormnum = '40';
strain = 'npr1';
marker = 'pharynx';
phase = 'fullMovie';
numTrajToPlot = 500;
applyMinTrajLength = false;
minTrajFrameNum = 45; % minimum number of frames required for the trajectory

if strcmp(marker,'pharynx')
    wormcats = {'loneWorm'}';
elseif strcmp(marker,'bodywall')
    wormcats = {'loneWorm','leaveCluster'};
end

% load files
load(['figures/turns/results/allTrajxcoords_' strain '_' wormnum '_' phase '_data' num2str(dataset) '_' marker '.mat'])
load(['figures/turns/results/allTrajycoords_' strain '_' wormnum '_' phase '_data' num2str(dataset) '_' marker '.mat'])

for wormcatCtr = 1:length(wormcats)
    % create figure to hold traj in
    eval([wormcats{wormcatCtr} 'TrajFig']) = figure; hold on
    % filter for minimum traj length, save up to 200 traj
    if applyMinTrajLength
        longTrajxcoords.(wormcats{wormcatCtr}) = cell(numTrajToPlot,1);
        longTrajycoords.(wormcats{wormcatCtr}) = cell(numTrajToPlot,1);
        longTrajxcoordsCtr.(wormcats{wormcatCtr}) = 1;
        longTrajycoordsCtr.(wormcats{wormcatCtr}) = 1;
        while longTrajxcoordsCtr.(wormcats{wormcatCtr})<numTrajToPlot
            for trajCtr = 1:length(allTrajxcoords.(wormcats{wormcatCtr}))
                currentTrajxcoords = allTrajxcoords.(wormcats{wormcatCtr}){trajCtr};
                if size(currentTrajxcoords,1) > minTrajFrameNum
                    longTrajxcoords.(wormcats{wormcatCtr}){longTrajxcoordsCtr.(wormcats{wormcatCtr})} = currentTrajxcoords;
                    longTrajycoords.(wormcats{wormcatCtr}){longTrajycoordsCtr.(wormcats{wormcatCtr})} = allTrajycoords.(wormcats{wormcatCtr}){trajCtr};
                    longTrajxcoordsCtr.(wormcats{wormcatCtr}) = longTrajxcoordsCtr.(wormcats{wormcatCtr})+1;
                    longTrajycoordsCtr.(wormcats{wormcatCtr}) = longTrajycoordsCtr.(wormcats{wormcatCtr})+1;
                end
            end
        end
    end
    % randomly sample a number of trajs
    if applyMinTrajLength
        trajInd = [1:length(longTrajxcoords.(wormcats{wormcatCtr}))];
    else
        trajInd = [1:length(allTrajxcoords.(wormcats{wormcatCtr}))];
    end
    sampleTrajInd = datasample(trajInd,numTrajToPlot-1,'Replace',false);
    trajLength.(wormcats{wormcatCtr}) = NaN(20,1);
    for trajCtr = 1:length(sampleTrajInd)
        % load xy coordinates
        if applyMinTrajLength
            xcoords = longTrajxcoords.(wormcats{wormcatCtr}){trajCtr};
            ycoords = longTrajycoords.(wormcats{wormcatCtr}){trajCtr};
        else
            xcoords = allTrajxcoords.(wormcats{wormcatCtr}){trajCtr};
            ycoords = allTrajycoords.(wormcats{wormcatCtr}){trajCtr};
        end
        % take centroid
        xcoords = mean(xcoords,2);
        ycoords = mean(ycoords,2);
        % set xy coordinates to start at 0
        xcoords = xcoords - xcoords(1);
        ycoords = ycoords - ycoords(1);
        % display length of traj
        trajLength.(wormcats{wormcatCtr})(trajCtr) = length(xcoords);
        % rotate the traj
        
        % plot traj
        plot(xcoords,ycoords)
    end
    set(0,'CurrentFigure', eval([wormcats{wormcatCtr} 'TrajFig']))
    title(['sample ' wormcats{wormcatCtr} ' trajectories, max ' num2str(max(trajLength.(wormcats{wormcatCtr}))) ' frames']);
    axis equal
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    xlabel('microns');
    ylabel('microns');
    if applyMinTrajLength
        xlim([-5000 5000])
        ylim([-5000 5000])
    end
end