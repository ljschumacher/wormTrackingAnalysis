% calculate various single worm statistics for different numbers of worms
% on plate

% issues / todo:
% - should distances be calculated only to other moving worms, or to any
% object (ie also worms that won't appear in the next frame)?

clear
close all

dataset = 1; % enter 1 or 2

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

pixelsize = 100/19.5; % 100 microns are 19.5 pixels

if dataset ==1
    strains = {'npr1','HA','N2'};
elseif dataset ==2
    strains = {'npr1','N2'};
end
wormnums = {'40','HD'};
maxBlobSize = 1e4;
intensityThresholds = containers.Map({'40','HD','1W'},{60, 40, 100});
clusterthreshold = 500;
plotDiagnostics = false;
if dataset ==1
    savesubfolder = 'green1/';
elseif dataset ==2
    savesubfolder = 'green2/';
end


for wormnum = wormnums
    clustFig = figure; hold on
    for strainCtr = 1:length(strains)
        %% load data
        if dataset ==1
            filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum{1} '_list.txt']);
        elseif dataset ==2
            filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum{1} '_g_list.txt']);
        end
        numFiles = length(filenames);
        nnDistances = cell(numFiles,1);
        nnnDistances = cell(numFiles,1);
        clustercount = cell(numFiles,1);
        for fileCtr = 1:numFiles
            filename = filenames{fileCtr};
            if exist(filename,'file')
                trajData = h5read(filename,'/trajectories_data');
                blobFeats = h5read(filename,'/blob_features');
                skelData = h5read(filename,'/skeleton');
                frameRate = h5readatt(filename,'/plate_worms','expected_fps');
                maxNumFrames = numel(unique(trajData.frame_number));
                numFrames = round(maxNumFrames/frameRate/10);
                framesAnalyzed = randperm(maxNumFrames,numFrames); % randomly sample frames without replacement
                %% filter worms
                trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                    intensityThresholds(wormnum{1}),maxBlobSize);
                %% calculate stats
                neighbr_distances = h5read(filename,'/neighbr_distances');
                nnDistances{fileCtr} = neighbr_distances(trajData.filtered,1);
                nnnDistances{fileCtr} = neighbr_distances(trajData.filtered,2);
                % count how many worms are within clusterthreshold
                num_close_neighbrs = h5read(filename,'/num_close_neighbrs');
                clustercount{fileCtr} = num_close_neighbrs(trajData.filtered);
            else
                warning(['Not all necessary tracking results present for ' filename ])
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
        figurename = ['nearestneighbrdistance_' wormnum{1} '_' strains{strainCtr}];
        exportfig(nnFig,['figures/neighbrDensity/' savesubfolder figurename '.eps'],exportOptions)
        system(['epstopdf figures/neighbrDensity/' savesubfolder figurename '.eps']);
        system(['rm figures/neighbrDensity/' savesubfolder figurename '.eps']);
    end
    %% format and export figures
    title(clustFig.Children,wormnum{1},'FontWeight','normal');
    set(clustFig,'PaperUnits','centimeters')
    clustFig.Children.XLim = [0 10];
    %     clustFig.Children.YLim = [0 0.5];
    xlabel(clustFig.Children,['# neighbours (n) within ' num2str(clusterthreshold) ' \mum'])
    ylabel(clustFig.Children,'cumulative P(N\geq 1)')
    legend(clustFig.Children,strains)
    figurename = ['clusthist_' wormnum{1}];
    exportfig(clustFig,['figures/neighbrDensity/' savesubfolder figurename '.eps'],exportOptions)
    system(['epstopdf figures/neighbrDensity/' savesubfolder figurename '.eps']);
    system(['rm figures/neighbrDensity/' savesubfolder figurename '.eps']);
end