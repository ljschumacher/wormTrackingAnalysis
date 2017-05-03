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

strains = {'npr1','N2','HA'};
wormnums = {'40','HD'};
intensityThresholds = [50, 40, 100];
maxBlobSize = 1e4;
plotDiagnostics = true;

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    distFig = figure; hold on
    for strainCtr = 1:length(strains)
        %% load data
        filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_list.txt']);
        numFiles = length(filenames);
        pairDistances = cell(numFiles,1);
        for fileCtr = 1:numFiles
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            frameRate = h5readatt(filename,'/plate_worms','expected_fps');
            maxNumFrames = numel(unique(trajData.frame_number));
            numFrames = round(maxNumFrames/frameRate/10);
            framesAnalyzed = randperm(maxNumFrames,numFrames); % randomly sample frames without replacement
            %% filter worms
            trajData.filtered = (blobFeats.area*pixelsize^2<=maxBlobSize)&...
                (blobFeats.intensity_mean>=intensityThresholds(numCtr));
            %% calculate stats
            pairDistances{fileCtr} = cell(numFrames,1);
            for frameCtr = 1:numFrames
                frame = framesAnalyzed(frameCtr);
                [x, y] = getWormPositions(trajData, frame, true);
                if numel(x)>1 % need at least two worms in frame to calculate distances
                    pairDistances{fileCtr}{frameCtr} = pdist([x y]).*pixelsize; % distance of every worm to every other
                end
            end
        end
        %% plot data
        % pool data from all frames for each file, then for all files
        pairDistances = cellfun(@(x) {horzcat(x{:})},pairDistances);
        histogram(horzcat(pairDistances{:}),'Normalization','Probability',...
            'DisplayStyle','stairs')
    end
    %% format and export figures
    title(distFig.Children,wormnum,'FontWeight','normal');
    set(distFig,'PaperUnits','centimeters')
    xlim([0 1.2e4])
    xlabel(distFig.Children,'pair distance (\mum)')
    ylabel(distFig.Children,'P')
    legend(strains)
    figurename = ['pairdistance_' wormnum];
    exportfig(distFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end