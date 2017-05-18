% plot small cluster worms with their neighbrs
close all
clear

% specify data sets
dataset = 2; % specify which dataset to run the script for. Enter either 1 or 2
strains = {'npr1','N2'};
wormnums = {'40','HD'};

% set parameters for filtering data
neighbrCutOff = 500; % distance in microns to consider a neighbr close
maxBlobSize = 1e4;
loneClusterRadius = 2000; % distance in microns to consider a cluster by itself
intensityThresholds = [60, 40];
pixelsize = 100/19.5; % 100 microns are 19.5 pixels

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',50,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        %% load data
        perduranceProbFig = figure;
        perduranceCountFig = figure;
        if dataset ==1
            filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_list.txt']);
        elseif dataset ==2
            filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_g_list.txt']);
        end
        numFiles = length(filenames);
        for fileCtr=1:numFiles
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            numCloseNeighbr = h5read(filename,'/num_close_neighbrs');
            neighbrDist = h5read(filename,'/neighbr_distances');
            %% filter data
            % filter green by blob size and intensity
            trajData.filtered = (blobFeats.area*pixelsize^2<=maxBlobSize)&...
                (blobFeats.intensity_mean>=intensityThresholds(numCtr));
            % filter green by small cluster status
            trajData.filtered = trajData.filtered&...
                ((numCloseNeighbr== 2 & neighbrDist(:,3)>=loneClusterRadius)...
                |(numCloseNeighbr== 3 & neighbrDist(:,4)>=(loneClusterRadius))...
                |(numCloseNeighbr== 4 & neighbrDist(:,5)>=(loneClusterRadius)));
            % generate distribution of small cluster frames
            smallClusterFrames = trajData.frame_number(trajData.filtered)';
            if isempty(smallClusterFrames) == false
                q = diff([0 diff([smallClusterFrames]) 0]==1);
                consFrames = find(q == -1) - find(q == 1) + 1; % list lengths of consecutive frames
                singleFrames = length(smallClusterFrames) - sum(consFrames(:)); % find number of single frames
                frameDist = [ones(1,singleFrames) consFrames]; % compile distribution of consecutive frames
                % plot data
                % percentage graph
                figure(perduranceProbFig);subplot(3,5,fileCtr)
                histogram(frameDist,'Normalization','Probability','DisplayStyle','bar','BinWidth',1)
                title (strrep(strrep(filename(end-32:end-17),'_',''),'/',''))
                xlabel('cluster perdurance (frames)')
                ylabel('probability')
                % count graph
                figure(perduranceCountFig);subplot(3,5,fileCtr)
                histogram(frameDist,'Normalization','count','DisplayStyle','bar','BinWidth',1)
                title (strrep(strrep(filename(end-32:end-17),'_',''),'/',''))
                xlabel('cluster perdurance (frames)')
                ylabel('count')
            end
        end
        %% format graphs and export
        figure(perduranceProbFig)
        set(perduranceProbFig,'Name',[strain ' ' wormnum ' '])
        if dataset ==1
            epsFileName = ['figures/smallClusterPerdurance/green1/pdf/smallClusterPerduranceProb_' strain '_' wormnum '.eps'];
            figFileName = ['figures/smallClusterPerdurance/green1/fig/smallClusterPerduranceProb_' strain '_' wormnum '.fig'];
        elseif dataset==2
            epsFileName = ['figures/smallClusterPerdurance/green2/pdf/smallClusterPerduranceProb_' strain '_' wormnum '.eps'];
            figFileName = ['figures/smallClusterPerdurance/green2/fig/smallClusterPerduranceProb_' strain '_' wormnum '.fig'];
        end
        savefig(figFileName)
        exportfig(perduranceProbFig,epsFileName,exportOptions)
        system(['epstopdf ' epsFileName]);
        system(['rm ' epsFileName]);
        %
        figure(perduranceCountFig);
        set(perduranceCountFig,'Name',[strain ' ' wormnum ' '])
        if dataset ==1
            epsFileName = ['figures/smallClusterPerdurance/green1/pdf/smallClusterPerduranceCount_' strain '_' wormnum '.eps'];
            figFileName = ['figures/smallClusterPerdurance/green1/fig/smallClusterPerduranceCount_' strain '_' wormnum '.fig'];
        elseif dataset ==2
            epsFileName = ['figures/smallClusterPerdurance/green2/pdf/smallClusterPerduranceCount_' strain '_' wormnum '.eps'];
            figFileName = ['figures/smallClusterPerdurance/green2/fig/smallClusterPerduranceCount_' strain '_' wormnum '.fig'];
        end
        savefig(figFileName)
        exportfig(perduranceCountFig,epsFileName,exportOptions)
        system(['epstopdf ' epsFileName]);
        system(['rm ' epsFileName]);
    end
end