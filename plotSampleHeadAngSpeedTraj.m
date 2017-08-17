function plotSampleHeadAngSpeedTraj(headAngSpeedRanges,wormcat,strain,numSampleTraj)

% import saved traj values
headAngSpeedSampleTraj = importdata('figures/turns/results/headAngSpeedSampleTraj_npr1_40.mat');
for rangeCtr = 1:size(headAngSpeedRanges,1)
    for wormcatCtr = 1:length(wormcat)
        samplePathFig = figure; hold on
        % remove empty cells
        headAngSpeedSampleTrajRange = squeeze(headAngSpeedSampleTraj.(wormcat{wormcatCtr})(:,:,rangeCtr));
        headAngSpeedSampleTrajRange =  headAngSpeedSampleTrajRange(~cellfun('isempty',headAngSpeedSampleTrajRange));
        headAngSpeedSampleTrajRange = reshape(headAngSpeedSampleTrajRange,[],2);
        % randomly sample 10 trajectories from the 500 saved ones
        trajSamples = randi(size(headAngSpeedSampleTrajRange,1),[numSampleTraj,1]);
        % loop through each saved trajectory xy coordinates
        for trajCtr = 1:length(trajSamples)
            % get xy coordinates for sample traj
            traj_xcoords = headAngSpeedSampleTraj.(wormcat{wormcatCtr}){trajSamples(trajCtr),1,rangeCtr};
            traj_ycoords = headAngSpeedSampleTraj.(wormcat{wormcatCtr}){trajSamples(trajCtr),2,rangeCtr};
            % take mean
            traj_xcoords = mean(traj_xcoords,2);
            traj_ycoords = mean(traj_ycoords,2);
            % set all trajectories to start at 0,0
            traj_xcoords = traj_xcoords - traj_xcoords(1);
            traj_ycoords = traj_ycoords - traj_ycoords(1);
            % plot
            set(0,'CurrentFigure',samplePathFig)
            plot(traj_xcoords,traj_ycoords)
            %waitforbuttonpress
        end
        % format and save plot
        title([strain '\_' wormcat{wormcatCtr} ' sampleTraj, '...
            num2str(headAngSpeedRanges(rangeCtr,1)) '-' num2str(headAngSpeedRanges(rangeCtr,2)) '°/s'],'FontWeight','normal')
        set(samplePathFig,'PaperUnits','centimeters')
        xlim([-100 100])
        ylim([-100 100])
        ax = gca;
        ax.XAxisLocation = 'origin'
        ax.YAxisLocation = 'origin'
        figurename = ['figures/turns/sampleTraj/' strain '_' wormcat{wormcatCtr} '_range' num2str(rangeCtr)];
        savefig(samplePathFig,[figurename '.fig'])
        load('exportOptions.mat')
        exportfig(samplePathFig,[figurename '.eps'],exportOptions)
        %system(['epstopdf ' figurename '.eps']);
        %system(['rm ' figurename '.eps']);
    end
end