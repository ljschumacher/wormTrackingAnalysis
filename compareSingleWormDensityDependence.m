% calculate various single worm statistics for different numbers of worms
% on plate

% issues / todo:
% - should distances be calculated only to other moving worms, or to any
% object (ie also worms that won't appear in the next frame)?

clear
close all

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

pixelsize = 100/19.5; % 100 microns are 19.5 pixels

strains = {'npr1','HA','N2'};
wormnums = {'40','HD','1W'};
intensityThresholds = containers.Map({'40','HD','1W'},{50, 40, 100});
maxBlobSize = 1e4;
plotDiagnostics = false;
numSamples = 8500;
maxSpeed = 1e3;
minNeighbrDist = 1500;

for strainCtr = 1:length(strains)
    speedFig = figure; hold on
    for wormnum = wormnums
        %% load data
        filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum{1} '_list.txt']);
        numFiles = length(filenames);
        speeds = cell(numFiles,1);
        loneWorms = cell(numFiles,1);
        for fileCtr = 1:numFiles % can be parfor
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            %% filter worms
            trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                    intensityThresholds(wormnum{1}),maxBlobSize);
            if ~strcmp(wormnum{1},'1W') % filter for lone worms
                min_neighbr_dist = h5read(filename,'/min_neighbr_dist');
                trajData.filtered = trajData.filtered&min_neighbr_dist>=minNeighbrDist;
            end
            %% select frames from filtered data only
            uniqueFrames = unique(trajData.frame_number(trajData.filtered));
            maxNumFrames = numel(uniqueFrames);
            numFrames = round(maxNumFrames/frameRate);
            framesAnalyzed = randsample(uniqueFrames,numFrames); % randomly sample frames without replacement
            %% calculate stats
            % beware that this ordering is frames first, not worms first
            % (as the saved files)
            speeds{fileCtr} = cell(numFrames,1);
            for frameCtr = 1:numFrames
                frame = framesAnalyzed(frameCtr);
                [~, ~, u, v] = calculateWormSpeeds(trajData, frame, true);
                speeds{fileCtr}{frameCtr} = sqrt(u.^2+v.^2)*pixelsize*frameRate; % speed of every worm in frame, in mu/s
            end
        end
        %% plot data
        % pool data from all frames for each file, then for all files
        speeds = cellfun(@(x) {vertcat(x{:})},speeds);
        speedscatenated = vertcat(speeds{:});
        % plot speed distribution
        speedscatenated = speedscatenated(speedscatenated<=maxSpeed);
        histogram(speedscatenated(randperm(numel(speedscatenated),numSamples)),...
                'DisplayStyle','stairs','Normalization','Probability','BinLimits',[0 maxSpeed])
    end
    %% format and export figures
    title(speedFig.Children,[strains{strainCtr} ', ' num2str(numSamples) ' samples'],'FontWeight','normal');
    set(speedFig,'PaperUnits','centimeters')
    xlabel(speedFig.Children,'speed (\mum/s)')
    ylabel(speedFig.Children,'P')
    legend(wormnums)
    figurename = ['figures/singleWorm/speeddistributions_' strains{strainCtr}];
    exportfig(speedFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end