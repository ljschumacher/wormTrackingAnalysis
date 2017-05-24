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
minSkelLength_r = 850;
maxSkelLength_r = 1500;
loneClusterRadius = 2000; % distance in microns to consider a cluster by itself
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
        filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_r_list.txt']);
        numFiles = length(filenames);
        frameDist = [];
        for fileCtr=1:numFiles
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton');
            numCloseNeighbr = h5read(filename,'/num_close_neighbrs');
            neighbrDist = h5read(filename,'/neighbr_distances');
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
            % generate distribution of small cluster frames
            smallClusterFrames = trajData.frame_number(trajData.filtered)';
            if isempty(smallClusterFrames) == false
                q = diff([0 diff([smallClusterFrames]) 0]==1);
                contFrames = find(q == -1) - find(q == 1) + 1; % list lengths of consecutive frames
                singleFrames = length(smallClusterFrames) - sum(contFrames(:)); % find number of single frames
                frameDist = [frameDist ones(1,singleFrames) contFrames]; % compile distribution of consecutive frames by adding each movie
            end
        end
        % plot cumulative survival
        [ecdfy,ecdfx] = ecdf(frameDist);
        plot(ecdfx,1-ecdfy)
        hold on
        legendMatrix{(numCtr-1)*2+(strainCtr)}= strcat(strain, '\_', wormnum);
    end
end
%% format graphs and export
xlabel('frames elapsed (at 9fps)')
xlim([0 40])
ylabel('remaining proportion')
legend(legendMatrix)
epsFileName = ['figures/smallClusterPerdurance/red2/pdf/smallClusterPerduranceSurvivalPooled.eps'];
figFileName = ['figures/smallClusterPerdurance/red2/fig/smallClusterPerduranceSurvivalPooled.fig'];
%savefig(figFileName)
%exportfig(cumSurvivalFig,epsFileName,exportOptions)
%system(['epstopdf ' epsFileName]);
%system(['rm ' epsFileName]);