% plot worm density over time so HD FOV density can be estimated for
% simulation purpose

clear
close all

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

%% set parameters
dataset = 1; % specify which dataset to run the script for. Enter either 1 or 2
if dataset ==1
    strains = {'N2','HA','npr1'};
elseif dataset ==2
    strains = {'N2','npr1'};
end
wormnums = {'1W','40','HD'};%{'1W','40','HD'};
if dataset == 1
    intensityThresholds = containers.Map({'40','HD','1W'},{50, 40, 100});
elseif dataset ==2
    intensityThresholds = containers.Map({'40','HD','1W'},{60, 40, 100});
end
maxBlobSize = 1e4;
pixelsize = 100/19.5; % 100 microns are 19.5 pixels

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        % load data
        if dataset == 1
            filenames = importdata(['../datalists/' strains{strainCtr} '_' wormnum '_list.txt']);
        elseif dataset == 2
            filenames = importdata(['../datalists/' strains{strainCtr} '_' wormnum '_g_list.txt']);
        end
        numFiles = length(filenames);
        wormDensityFig = figure; hold on
        for fileCtr=1:numFiles
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            % filter green by blob size and intensity
            trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                intensityThresholds(wormnum),maxBlobSize);
            % calculate objects per frame
            objInFrame = trajData.frame_number(trajData.filtered);
            % plot histogram
            set(0,'CurrentFigure',wormDensityFig)
            plotColor = colorcube(numFiles);
            histogram(objInFrame,'BinWidth',45,'DisplayStyle','stairs','EdgeColor',plotColor(fileCtr,:),'Normalization','countdensity')
        end
        % format figure and export
        xlabel('frame number','FontSize',20)
        ylabel('number of tracked objects','FontSize',20)
        set(gca,'FontSize',15)
        title([strains{strainCtr} ' ' wormnums{numCtr} ' data ' num2str(dataset)],'FontWeight','normal');
        figurename = ['../figures/wormDensity/wormDensity_' strains{strainCtr} '_' wormnums{numCtr} '_data' num2str(dataset)];
        exportfig(wormDensityFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
    end
end