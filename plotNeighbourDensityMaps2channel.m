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
maxBlobSize_r = 2.5e5;
minSkelLength = 850;
maxSkelLength = 1500;
clusterthreshold = 500;
plotDiagnostics = false;

for wormnum = wormnums
    clustFig = figure; hold on
    for strainCtr = 1:length(strains)
        %% load data
        filenames_r = importdata(['datalists/' strains{strainCtr} '_' wormnum{1} '_r_list.txt']);
        numFiles = length(filenames_r);
        assert(length(filenames_r)==numFiles,'Number of files for two channels do not match.')
        nnDistances = cell(numFiles,1);
        nnnDistances = cell(numFiles,1);
        clustercount = cell(numFiles,1);
        for fileCtr = 1:numFiles
            filename_r = filenames_r{fileCtr};
            if exist(filename_r,'file')
                trajData_r = h5read(filename_r,'/trajectories_data');
                blobFeats_r = h5read(filename_r,'/blob_features');
                skelData_r = h5read(filename_r,'/skeleton');
                frameRate = h5readatt(filename_r,'/plate_worms','expected_fps');
                maxNumFrames = numel(unique(trajData_r.frame_number));
                numFrames = round(maxNumFrames/frameRate/10);
                framesAnalyzed = randperm(maxNumFrames,numFrames); % randomly sample frames without replacement
                %% filter worms
                if contains(filename_r,'55')||contains(filename_r,'54')
                    intensityThreshold_r = 80;
                else
                    intensityThreshold_r = 40;
                end
                trajData_r.filtered = filterIntensityAndSize(blobFeats_r,pixelsize,...
                    intensityThreshold_r,maxBlobSize_r)&...
                    logical(trajData_r.is_good_skel)&...
                filterSkelLength(skelData_r,pixelsize,minSkelLength,maxSkelLength);
                %% calculate stats
                neighbr_distances = h5read(filename_r,'/neighbr_distances');
                nnDistances{fileCtr} = neighbr_distances(trajData_r.filtered,1);
                nnnDistances{fileCtr} = neighbr_distances(trajData_r.filtered,2);
                % count how many worms are within clusterthreshold
                num_close_neighbrs = h5read(filename_r,'/num_close_neighbrs');
                clustercount{fileCtr} = num_close_neighbrs(trajData_r.filtered);
            else
                warning(['Not all necessary tracking results present for ' filename_r ])
                nnDistances{fileCtr} = [];
                nnnDistances{fileCtr} = [];
            end
        end
        %% plot data
        % pool data from all frames for each file, then for all files
        ecdf(clustFig.Children,vertcat(clustercount{:}),'function','survivor')
        % nearest vs next-nearest neighbours
        nnDistances = cellfun(@(x) {vertcat(x)},nnDistances);
        nnnDistances = cellfun(@(x) {vertcat(x)},nnnDistances);
        nnFig = figure;
        histogram2(vertcat(nnDistances{:}),vertcat(nnnDistances{:}),...
            'BinWidth',clusterthreshold*[1 1]/2,'Normalization','Probability','DisplayStyle','tile','EdgeColor','none')
        xlim([0 4000]), ylim([0 4000])
        title(nnFig.Children,[wormnum{1} ' ' strains{strainCtr}],'FontWeight','normal');
        set(nnFig,'PaperUnits','centimeters')
        xlabel(nnFig.Children,'nn distance (\mum)')
        ylabel(nnFig.Children,'nnn distance (\mum)')
        figurename = ['nearestneighbrdistance_rg_' wormnum{1} '_' strains{strainCtr}];
        exportfig(nnFig,['figures/' figurename '.eps'],exportOptions)
        system(['epstopdf figures/' figurename '.eps']);
        system(['rm figures/' figurename '.eps']);
    end
    %% format and export figures 
    title(clustFig.Children,wormnum{1},'FontWeight','normal');
    set(clustFig,'PaperUnits','centimeters')
    clustFig.Children.XLim = [0 10];
%     clustFig.Children.YLim = [0 0.5];
    xlabel(clustFig.Children,['# neighbours (n) within ' num2str(clusterthreshold) ' \mum'])
    ylabel(clustFig.Children,'cumulative P(N\geq 1)')
    legend(clustFig.Children,strains)
    figurename = ['clusthist_rg_' wormnum{1}];
    exportfig(clustFig,['figures/' figurename '.eps'],exportOptions)
    system(['epstopdf figures/' figurename '.eps']);
    system(['rm figures/' figurename '.eps']);
end