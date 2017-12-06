% plot small cluster worms with their neighbrs
close all
clear

% specify data sets
dataset = 1; % specify which dataset to run the script for. Enter either 1 or 2
if dataset == 1
    strains = {'npr1','HA','N2'};
elseif dataset ==2
    strains = {'npr1','N2'};
end
wormnums = {'40','HD'};

% set parameters for filtering data
neighbrCutOff = 500; % distance in microns to consider a neighbr close
maxBlobSize = 1e4;
minNeighbrDist = 2000; % distance in microns to consider a cluster by itself
intensityThresholds = containers.Map({'40','HD','1W'},{60, 40, 100});
pixelsize = 100/19.5; % 100 microns are 19.5 pixels

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',30,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        %% load data
        if dataset ==1
            filenames = importdata(['../datalists/' strains{strainCtr} '_' wormnum '_list.txt']);
        elseif dataset ==2
            filenames = importdata(['../datalists/' strains{strainCtr} '_' wormnum '_g_list.txt']);
        end
        numFiles = length(filenames);
        if numFiles >12
            numFiles = 12; % only plot max 12 individual movies
        end
        speedDistFig = figure;
        for fileCtr=1:numFiles
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            min_neighbr_dist = h5read(filename,'/min_neighbr_dist');
            num_close_neighbrs = h5read(filename,'/num_close_neighbrs');
            neighbr_dist = h5read(filename,'/neighbr_distances');
            %% filter data
            % filter green by blob size and intensity
            trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                intensityThresholds(wormnum),maxBlobSize);
            % define cluster status
            all = trajData.filtered;
            inCluster = trajData.filtered& num_close_neighbrs>=3;
            smallCluster = trajData.filtered&...
                ((num_close_neighbrs == 2 & neighbr_dist(:,3)>=(minNeighbrDist))...
                |(num_close_neighbrs == 3 & neighbr_dist(:,4)>=(minNeighbrDist))...
                |(num_close_neighbrs == 4 & neighbr_dist(:,5)>=(minNeighbrDist)));
            loneWorms = min_neighbr_dist>=minNeighbrDist;
            % generate speed distribution for each cluster class
            allSpeeds = blobFeats.signed_speed(all);
            inClusterSpeeds = blobFeats.signed_speed(inCluster);
            smallClusterSpeeds = blobFeats.signed_speed(smallCluster);
            loneWormSpeeds = blobFeats.signed_speed(loneWorms);
            % trim speed distribution down to between -15 and 15
            allSpeeds = allSpeeds(allSpeeds<=15&allSpeeds>=-15);
            inClusterSpeeds = inClusterSpeeds(inClusterSpeeds<=15&inClusterSpeeds>=-15);
            smallClusterSpeeds = smallClusterSpeeds(smallClusterSpeeds<=15&smallClusterSpeeds>=-15);
            loneWormSpeeds = loneWormSpeeds(loneWormSpeeds<=15&loneWormSpeeds>=-15);
            % plot speed distribution
            individualplot = subplot(3,4,fileCtr); hold on
            histogram(allSpeeds,'Normalization','probability','DisplayStyle','stairs','BinWidth',0.5)
            histogram(inClusterSpeeds,'Normalization','probability','DisplayStyle','stairs','BinWidth',0.5)
            histogram(smallClusterSpeeds,'Normalization','probability','DisplayStyle','stairs','BinWidth',0.5)
            histogram(loneWormSpeeds,'Normalization','probability','DisplayStyle','stairs','BinWidth',0.5)
            xlabel('speed')
            ylabel('probability')
            individualplot.XLim = [-15 15];
            individualplot.YLim = [0 0.25];
            title(individualplot,strrep(strrep(filename(end-32:end-17),'_',''),'/',''))
            legend('all','in cluster','small cluster','lone worms','Location','Northwest')
        end
        % save figure
        if dataset ==1
            epsFileName = ['../figures/speeds/g1_' strains{strainCtr} '_' wormnum '.eps'];
            figFileName = ['../figures/speeds/g1_' strains{strainCtr} '_' wormnum '.fig'];
        elseif dataset ==2
            epsFileName = ['../figures/speeds/g2_' strains{strainCtr} '_' wormnum '.eps'];
            figFileName = ['../figures/speeds/g2_' strains{strainCtr} '_' wormnum '.fig'];
        end
        savefig(figFileName)
        exportfig(speedDistFig,epsFileName,exportOptions)
        system(['epstopdf ' epsFileName]);
        system(['rm ' epsFileName]);
    end
end