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

strains = {'npr1','N2'};
wormnums = {'40','HD'};
intensityThresholds_g = [60, 40];
maxBlobSize_g = 1e4;
maxBlobSize_r = 2.5e5;
clusterthreshold = 1000;
plotDiagnostics = false;

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    distFig = figure; hold on
    clustFig = figure; hold on
    for strainCtr = 1:length(strains)
        %% load data
        filenames_g = importdata(['datalists/' strains{strainCtr} '_' wormnum '_g_list.txt']);
        filenames_r = importdata(['datalists/' strains{strainCtr} '_' wormnum '_r_list.txt']);
        numFiles = length(filenames_g);
        assert(length(filenames_r)==numFiles,'Number of files for two channels do not match.')
        pairDistances = cell(numFiles,1);
        nnDistances = cell(numFiles,1);
        nnnDistances = cell(numFiles,1);
        clustercount = cell(numFiles,1);
        for fileCtr = 1:numFiles
            filename_g = filenames_g{fileCtr};
            filename_r = filenames_r{fileCtr};
            if exist(filename_r,'file')&&exist(filename_g,'file')
                trajData_g = h5read(filename_g,'/trajectories_data');
                trajData_r = h5read(filename_r,'/trajectories_data');
                blobFeats_g = h5read(filename_g,'/blob_features');
                blobFeats_r = h5read(filename_r,'/blob_features');
                frameRate = h5readatt(filename_g,'/plate_worms','expected_fps');
                maxNumFrames = numel(unique(trajData_g.frame_number));
                numFrames = round(maxNumFrames/frameRate/10);
                framesAnalyzed = randperm(maxNumFrames,numFrames); % randomly sample frames without replacement
                %% filter worms
                if contains(filename_r,'55')||contains(filename_r,'54')
                    intensityThreshold_r = 80;
                else
                    intensityThreshold_r = 40;
                end
                trajData_g.filtered = filterIntensityAndSize(blobFeats_g,pixelsize,...
                    intensityThresholds_g(numCtr),maxBlobSize_g);
                trajData_r.filtered = filterIntensityAndSize(blobFeats_r,pixelsize,...
                    intensityThreshold_r,maxBlobSize_r)&...
                    logical(trajData_r.is_good_skel);
                %% calculate stats
                pairDistances{fileCtr} = cell(numFrames,1);
                nnDistances{fileCtr} = cell(numFrames,1);
                nnnDistances{fileCtr} = cell(numFrames,1);
                for frameCtr = 1:numFrames
                    frame = framesAnalyzed(frameCtr);
                    [x_g, y_g] = getWormPositions(trajData_g, frame);
                    [x_r, y_r] = getWormPositions(trajData_r, frame);
                    if numel(x_g)>=1&&numel(x_r)>=1 % need at least two worms in frame to calculate distances
                        redToGreenDistances = sort(pdist2([x_r y_r],[x_g y_g]).*pixelsize,2); % distance of every red worm to every green
                        pairDistances{fileCtr}{frameCtr} = redToGreenDistances(:);
                        nnDistances{fileCtr}{frameCtr} = redToGreenDistances(:,1); % nearest neighbours
                        if size(redToGreenDistances,2)>1
                            nnnDistances{fileCtr}{frameCtr} = redToGreenDistances(:,2); % next-nearest neighbours
                        else
                            nnnDistances{fileCtr}{frameCtr} = NaN(size(x_r));
                        end
                    else
                        pairDistances{fileCtr}{frameCtr} = NaN;
                        nnDistances{fileCtr}{frameCtr} = NaN;
                        nnnDistances{fileCtr}{frameCtr} = NaN;
                    end
                end
                % count how many worms are within clusterthreshold
                clustercount{fileCtr} = cellfun(@(x) nnz(x<=clusterthreshold), pairDistances{fileCtr});
            else
                warning(['Not all necessary tracking results present for ' filename_r ])
                pairDistances{fileCtr} = {};
                nnDistances{fileCtr} = {};
                nnnDistances{fileCtr} = {};
            end
        end
        %% plot data
        % pool data from all frames for each file, then for all files
        pairDistances = cellfun(@(x) {vertcat(x{:})},pairDistances);
        histogram(distFig.Children,vertcat(pairDistances{:}),'Normalization','Probability',...
            'DisplayStyle','stairs')
        histogram(clustFig.Children,vertcat(clustercount{:}),'Normalization','Probability',...
            'DisplayStyle','stairs')
        % nearest vs next-nearest neighbours
        nnDistances = cellfun(@(x) {vertcat(x{:})},nnDistances);
        nnnDistances = cellfun(@(x) {vertcat(x{:})},nnnDistances);
        nnFig = figure;
        histogram2(vertcat(nnDistances{:}),vertcat(nnnDistances{:}),...
            'Normalization','Probability','DisplayStyle','tile','EdgeColor','none')
        xlim([0 4000]), ylim([0 4000])
        title(nnFig.Children,[wormnum ' ' strains{strainCtr}],'FontWeight','normal');
        set(nnFig,'PaperUnits','centimeters')
        xlabel(nnFig.Children,'nn distance (\mum)')
        ylabel(nnFig.Children,'nnn distance (\mum)')
        figurename = ['nearestneighbourdistance_rg_' wormnum '_' strains{strainCtr}];
        exportfig(nnFig,['figures/' figurename '.eps'],exportOptions)
        system(['epstopdf figures/' figurename '.eps']);
        system(['rm figures/' figurename '.eps']);
    end
    %% format and export figures
    title(distFig.Children,wormnum,'FontWeight','normal');
    set(distFig,'PaperUnits','centimeters')
    distFig.Children.XLim = [0 1.2e4];
    xlabel(distFig.Children,'r-g pair distance (\mum)')
    ylabel(distFig.Children,'P')
    legend(distFig.Children,strains)
    figurename = ['pairdistance_rg_' wormnum];
    exportfig(distFig,['figures/' figurename '.eps'],exportOptions)
    system(['epstopdf figures/' figurename '.eps']);
    system(['rm figures/' figurename '.eps']);
    %
    title(clustFig.Children,wormnum,'FontWeight','normal');
    set(clustFig,'PaperUnits','centimeters')
    clustFig.Children.XLim = [0 40];
    clustFig.Children.YLim = [0 0.25];
    xlabel(clustFig.Children,['# neighbours within ' num2str(clusterthreshold) ' \mum'])
    ylabel(clustFig.Children,'P')
    legend(clustFig.Children,strains)
    figurename = ['clusthist_rg_' wormnum];
    exportfig(clustFig,['figures/' figurename '.eps'],exportOptions)
    system(['epstopdf figures/' figurename '.eps']);
    system(['rm figures/' figurename '.eps']);
end