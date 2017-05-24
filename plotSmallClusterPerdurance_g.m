% plot small cluster worms with their neighbrs
close all
clear

% specify data sets
dataset = 1; % specify which dataset to run the script for. Enter either 1 or 2
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

cumSurvivalFig = figure;
legendMatrix=cell(length(strains)*length(wormnums),1);

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        %% load data
        if dataset ==1
            filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_list.txt']);
        elseif dataset ==2
            filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_g_list.txt']);
        end
        numFiles = length(filenames);
        frameDist = [];
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
                frameDist = [frameDist ones(1,singleFrames) consFrames]; % compile distribution of consecutive frames by adding each movie
            end
        end
        % plot cumulative survival
        [ecdfy,ecdfx] = ecdf(frameDist);
        plot(ecdfx,1-ecdfy)
        %ecdf(frameDist,'function','survivor','alpha',0.01,'bounds','on')
        hold on
        legendMatrix{(numCtr-1)*2+(strainCtr)}= strcat(strain, '\_', wormnum);
    end
end
%% format graphs and export
xlabel('frames elapsed (at 9fps)')
xlim([0 40])
ylabel('remaining proportion')
legend(legendMatrix)
if dataset ==1
    epsFileName = ['figures/smallClusterPerdurance/green1/pdf/smallClusterPerduranceSurvivalPooled.eps'];
    figFileName = ['figures/smallClusterPerdurance/green1/fig/smallClusterPerduranceSurvivalPooled.fig'];
elseif dataset ==2
    epsFileName = ['figures/smallClusterPerdurance/green2/pdf/smallClusterPerduranceSurvivalPooled.eps'];
    figFileName = ['figures/smallClusterPerdurance/green2/fig/smallClusterPerduranceSurvivalPooled.fig'];
end
savefig(figFileName)
exportfig(cumSurvivalFig,epsFileName,exportOptions)
system(['epstopdf ' epsFileName]);
system(['rm ' epsFileName]);