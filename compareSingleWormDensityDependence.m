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

strains = {'HA','npr1','N2'};
wormnums = {'HD','40','1W'};
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
        filenames = importdata([strains{strainCtr} '_' wormnum '_list.txt']);
        numFiles = length(filenames);
        speeds = cell(numFiles,1);
        mindist = cell(numFiles,1);
        parfor fileCtr = 1:numFiles % can be parfor
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            frameRate = h5readatt(filename,'/plate_worms','expected_fps');
            maxNumFrames = numel(unique(trajData.frame_number));
            numFrames = round(maxNumFrames/frameRate);
            framesAnalyzed = randperm(maxNumFrames,numFrames); % randomly sample frames without replacement
            %% filter worms
            if plotDiagnostics
                plotIntensitySizeFilter(blobFeats,pixelsize,...
                    intensityThresholds(numCtr),maxBlobSize,...
                    [wormnum ' ' strains{strainCtr} ' ' strrep(filename(end-38:end-23),'/','')])
            end
            trajData.filtered = (blobFeats.area*pixelsize^2<=maxBlobSize)&...
                (blobFeats.intensity_mean>=intensityThresholds(numCtr));
            %% calculate stats
            speeds{fileCtr} = cell(numFrames,1);
            mindist{fileCtr} = cell(numFrames,1);
            for frameCtr = 1:numFrames
                frame = framesAnalyzed(frameCtr);
                [x, y, u, v] = calculateWormSpeeds(trajData, frame);
                speeds{fileCtr}{frameCtr} = sqrt(u.^2+v.^2)*pixelsize*frameRate; % speed of every worm in frame, in mu/s
                if numel(x)>1 % need at least two worms in frame to calculate distances
                    D = squareform(pdist([x y]).*pixelsize); % distance of every worm to every other
                    mindist{fileCtr}{frameCtr} = min(D + max(max(D))*eye(size(D)));
                    if (numel(speeds{fileCtr}{frameCtr})~=numel(mindist{fileCtr}{frameCtr}))
                        error(['Inconsistent number of variables in frame ' num2str(frame) ' of ' filename ])
                    end
                elseif numel(x)==1
                    mindist{fileCtr}{frameCtr} = NaN;
                end
            end
        end
        %% plot data
        % pool data from all frames for each file, then for all files
        speeds = cellfun(@(x) {vertcat(x{:})},speeds);
        mindist = cellfun(@(x) {horzcat(x{:})},mindist);
        speedscatenated = vertcat(speeds{:});
        if ~strcmp(wormnum,'1W')
            % if it is multiworm data, we need to filter for isolated worms
            loneWorms = horzcat(mindist{:})>=2500;  
            speedscatenated = speedscatenated(loneWorms);
        end
        % filter for speeds
        speedscatenated = speedscatenated(speedscatenated<=maxSpeed);
        histogram(speedscatenated(randperm(numel(speedscatenated),numSamples)),...
                'DisplayStyle','stairs','Normalization','Probability','BinLimits',[0 maxSpeed])
        %% also fix bins... or pool data from multiple recordings
    end
    %% format and export figures
    title(speedFig.Children,[strains{strainCtr} ', ' num2str(numSamples) ' samples'],'FontWeight','normal');
    set(speedFig,'PaperUnits','centimeters')
    xlabel(speedFig.Children,'speed (\mum/s)')
    ylabel(speedFig.Children,'P')
    legend(wormnums)
    figurename = ['singleWorm/' strains{strainCtr} '_speeddistributions'];
    exportfig(speedFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end