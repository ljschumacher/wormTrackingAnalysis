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
intensityThresholds = [50, 40, 100];
maxBlobSize = 1e4;
plotDiagnostics = false;
numSamples = 8500;
maxSpeed = 1e3;

for strainCtr = 1:length(strains)
    speedFig = figure; hold on
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        %% load data
        filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_list.txt']);
        numFiles = length(filenames);
        speeds = cell(numFiles,1);
        loneWorms = cell(numFiles,1);
        for fileCtr = 1:numFiles % can be parfor
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            frameRate = h5readatt(filename,'/plate_worms','expected_fps');
            maxNumFrames = numel(unique(trajData.frame_number));
            numFrames = round(maxNumFrames/frameRate);
            framesAnalyzed = randperm(maxNumFrames,numFrames); % randomly sample frames without replacement
            %% filter worms
            trajData.filtered = (blobFeats.area*pixelsize^2<=maxBlobSize)&...
                (blobFeats.intensity_mean>=intensityThresholds(numCtr));
            %% calculate stats
            % beware that this ordering is frames first, not worms first
            % (as the saved files)
            speeds{fileCtr} = cell(numFrames,1);
            loneWorms{fileCtr} = cell(numFrames,1);
            for frameCtr = 1:numFrames
                frame = framesAnalyzed(frameCtr);
                [~, ~, u, v] = calculateWormSpeeds(trajData, frame);
                speeds{fileCtr}{frameCtr} = sqrt(u.^2+v.^2)*pixelsize*frameRate; % speed of every worm in frame, in mu/s
% %                 [~, loneWorms{fileCtr}{frameCtr}, ~] = ...
% %                     getWormClusterStatus(trajData, frame, pixelsize, 2500,500,2);
                if (numel(speeds{fileCtr}{frameCtr})~=numel(loneWorms{fileCtr}{frameCtr}))
                    error(['Inconsistent number of variables in frame ' num2str(frame) ' of ' filename ])
                end
            end
        end
        %% plot data
        % pool data from all frames for each file, then for all files
        speeds = cellfun(@(x) {vertcat(x{:})},speeds);
        loneWorms = cellfun(@(x) {horzcat(x{:})},loneWorms);
        speedscatenated = vertcat(speeds{:});
        if ~strcmp(wormnum,'1W')
            % if it is multiworm data, we need to filter for isolated worms
            loneWorms = horzcat(loneWorms{:});  
            speedscatenated = speedscatenated(loneWorms);
        end
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
    figurename = ['figures/singleWorm/' strains{strainCtr} '_speeddistributions'];
    exportfig(speedFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end