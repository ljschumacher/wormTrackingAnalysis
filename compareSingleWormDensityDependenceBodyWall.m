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
plotDiagnostics = false;
numSamples = 8500;

for strainCtr = 1:length(strains)
    speedFig = figure; hold on
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        %% load data
        filenames = importdata([strains{strainCtr} '_' wormnum '_r_list.txt']);
        numFiles = length(filenames);
        speeds = cell(numFiles,1);
        mindist = cell(numFiles,1);
        for fileCtr = 1:numFiles % can be parfor
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            featData = h5read(strrep(filename,'skeletons','features'),'/features_timeseries');
            frameRate = h5readatt(filename,'/plate_worms','expected_fps');
            maxNumFrames = numel(unique(trajData.frame_number));
            numFrames = round(maxNumFrames/frameRate);
            framesAnalyzed = randperm(maxNumFrames,numFrames); % randomly sample frames without replacement
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
            % find for which worms features have been calculated
            featData.filtered = ismember(int32(features.worm_index),...
                unique(trajData.worm_index_joined(trajData.filtered)));
% % % %             % filter features for whether worms are lone or in cluster?
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