% plot small cluster worms with their neighbrs
close all
clear

% specify data sets
strains = {'npr1','N2'};
wormnums = {'40','HD'};
numFramesSampled = 10; % how many frames to randomly sample per file

% set parameters for filtering data
neighbrCutOff = 500; % distance in microns to consider a neighbr close
maxBlobSize_r = 2.5e5;
maxBlobSize_g = 1e4;
minSkelLength_r = 850;
maxSkelLength_r = 1500;
loneClusterRadius = 2000; % distance in microns to consider a cluster by itself
intensityThresholds_g = [60, 40, NaN];
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
        filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_r_list.txt']);
        numFiles = length(filenames);
        filenames_g = importdata(['datalists/' strains{strainCtr} '_' wormnum '_g_list.txt']);
        perduranceProbFig = figure;
        perduranceCountFig = figure;
        for fileCtr=1:numFiles
            filename = filenames{fileCtr};
            filename_g = filenames_g{fileCtr};
            if exist(filename,'file')&&exist(filename_g,'file')
                trajData = h5read(filename,'/trajectories_data');
                blobFeats = h5read(filename,'/blob_features');
                skelData = h5read(filename,'/skeleton');
                numCloseNeighbr = h5read(filename,'/num_close_neighbrs');
                neighbrDist = h5read(filename,'/neighbr_distances');
                trajData_g = h5read(filename_g,'/trajectories_data');
                blobFeats_g = h5read(filename_g,'/blob_features');
                frameRate = h5readatt(filename,'/plate_worms','expected_fps');
                %% filter data
                % filter red by blob size and intensity
                if contains(filename,'55')||contains(filename,'54')
                    intensityThreshold = 80;
                else
                    intensityThreshold = 40;
                end
                trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                    intensityThreshold,maxBlobSize_r);
                % filter red by skeleton length
                trajData.filtered = trajData.filtered&logical(trajData.is_good_skel)...
                    &filterSkelLength(skelData,pixelsize,minSkelLength_r,maxSkelLength_r);
                % filter red by small cluster status
                trajData.filtered = trajData.filtered&...
                    ((numCloseNeighbr== 2 & neighbrDist(:,3)>=(loneClusterRadius))...
                    |(numCloseNeighbr== 3 & neighbrDist(:,4)>=(loneClusterRadius))...
                    |(numCloseNeighbr== 4 & neighbrDist(:,5)>=(loneClusterRadius)));
                % filter green channel by blob size and intensity
                trajData_g.filtered = (blobFeats_g.area*pixelsize^2<=maxBlobSize_g)&...
                    (blobFeats_g.intensity_mean>=intensityThresholds_g(numCtr));
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
                    ymax = length(find(frameDist == mode(frameDist)));
                    yticks(0:1:ymax)
                end
            end
        end
        %% format graphs and export
        figure(perduranceProbFig)
        set(perduranceProbFig,'Name',[strain ' ' wormnum ' '])
        epsFileName = ['figures/smallClusterPerdurance/red2/pdf/smallClusterPerduranceProb_' strain '_' wormnum '.eps'];
        figFileName = ['figures/smallClusterPerdurance/red2/pdf/smallClusterPerduranceProb_' strain '_' wormnum '.fig'];
        savefig(figFileName)
        exportfig(perduranceProbFig,epsFileName,exportOptions)
        system(['epstopdf ' epsFileName]);
        system(['rm ' epsFileName]);
        %
        figure(perduranceCountFig);
        set(perduranceCountFig,'Name',[strain ' ' wormnum ' '])
        epsFileName = ['figures/smallClusterPerdurance/red2/pdf/smallClusterPerduranceCount_' strain '_' wormnum '.eps'];
        figFileName = ['figures/smallClusterPerdurance/red2/pdf/smallClusterPerduranceCount_' strain '_' wormnum '.fig'];
        savefig(figFileName)
        exportfig(perduranceCountFig,epsFileName,exportOptions)
        system(['epstopdf ' epsFileName]);
        system(['rm ' epsFileName]);
    end
end