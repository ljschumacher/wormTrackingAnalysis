% plot in cluster worms
close all
clear

% specify data sets
strains = {'npr1','N2'};
wormnums = {'40','HD'};
numFramesSampled = 10; % how many frames to randomly sample per file

% set parameters for filtering data
neighbrCutOff = 500; % distance in microns to consider a neighbr close
maxBlobSize = 2.5e5;
maxBlobSize_g = 1e4;
minSkelLength = 850;
maxSkelLength = 1500;
loneClusterRadius = 2000; % distance in microns to consider a cluster by itself
intensityThresholds_g = [60, 40, NaN];
pixelsize = 100/19.5; % 100 microns are 19.5 pixels

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        %% load data
        filenames = importdata([strains{strainCtr} '_' wormnum '_r_list.txt']);
        numFiles = length(filenames);
        filenames_g = importdata([strains{strainCtr} '_' wormnum '_g_list.txt']);
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
                %% filter data
                % filter by blob size and intensity
                if contains(filename,'55')||contains(filename,'54')
                    intensityThreshold = 80;
                else
                    intensityThreshold = 40;
                end
                trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                    intensityThreshold,maxBlobSize);
                % filter by skeleton length
                trajData.filtered = trajData.filtered&logical(trajData.is_good_skel)...
                    &filterSkelLength(skelData,pixelsize,minSkelLength,maxSkelLength);
                % filter by small cluster status
                trajData.filtered = trajData.filtered&...
                    ((numCloseNeighbr== 2 & neighbrDist(:,3)>=(loneClusterRadius))...
                    |(numCloseNeighbr== 3 & neighbrDist(:,4)>=(loneClusterRadius))...
                    |(numCloseNeighbr== 4 & neighbrDist(:,5)>=(loneClusterRadius)));
                % filter green channel by blob size and intensity
                trajData_g.filtered = (blobFeats_g.area*pixelsize^2<=maxBlobSize_g)&...
                    (blobFeats_g.intensity_mean>=intensityThresholds_g(numCtr));
                % plot sample data
                if length(unique(trajData.frame_number(trajData.filtered)))<numFramesSampled
                    warning(['Not enough frames to plot for ' filename ])
                elseif unique(trajData.frame_number(trajData.filtered))>=numFramesSampled
                    smallClusterWormsFig = figure;
                    framesAnalyzed = randsample(unique(trajData.frame_number(trajData.filtered)),numFramesSampled);
                    for frameCtr = 1:numFramesSampled
                        frame=framesAnalyzed(frameCtr);
                        subplot(2,5,frameCtr)
                        % plot red worm
                        frameIdcs_worm = find(trajData.frame_number==frame&trajData.filtered);
                        if nnz(frameIdcs_worm)>1 % plot only one red worm if multiple present
                            frameIdcs_worm = randsample(frameIdcs_worm,1);
                        end
                        worm_xcoords = squeeze(skelData(1,:,frameIdcs_worm));
                        worm_ycoords = squeeze(skelData(2,:,frameIdcs_worm));
                        plot(worm_xcoords,worm_ycoords,'LineWidth',4)
                        hold on
                        % plot green worms
                        frameLogIdcs_pharynx = trajData_g.frame_number==frame&trajData_g.filtered;
                        plot(trajData_g.coord_x(frameLogIdcs_pharynx),trajData_g.coord_y(frameLogIdcs_pharynx),...
                            'ro','MarkerSize',4,'MarkerFaceColor','r')
                        axis equal
                        % plot circle of radius neighbrCutOff around each worm
                        viscircles([trajData.coord_x(frameIdcs_worm) trajData.coord_y(frameIdcs_worm)],...
                            neighbrCutOff/pixelsize,'LineStyle','--','Color',0.5*[1 1 1],'EnhanceVisibility',false);
                        % plot circle of radius loneClusterRadius around each worm
                        viscircles([trajData.coord_x(frameIdcs_worm) trajData.coord_y(frameIdcs_worm)],...
                            loneClusterRadius/pixelsize,'LineStyle','--','Color',0.5*[1 1 1],'EnhanceVisibility',false);
                        % plot circle of radius 2500 around each worm
                        viscircles([trajData.coord_x(frameIdcs_worm) trajData.coord_y(frameIdcs_worm)],...
                            2500/pixelsize,'LineStyle','--','Color',0.5*[1 1 1],'EnhanceVisibility',false);
                        %%% make separate plots for worms in clusters / lone
                        ax = gca;
                        xlim([-2500 2500]/pixelsize + trajData.coord_x(frameIdcs_worm))
                        ylim([-2500 2500]/pixelsize + trajData.coord_y(frameIdcs_worm))
                        set(ax,'visible','off')
                        movieName = strrep(strrep(filename(end-32:end-17),'_',''),'/','');
                        title = strcat(strain,'_',wormnum,'_',movieName);
                        ax.Position = ax.Position.*[1 1 1.2 1.2]; % reduce whitespace btw subplots
                    end
                    %% export figure
                    figName = strcat('smallCluster_',num2str(loneClusterRadius), '_',...
                        strain,'_',wormnum,'_',movieName,'.fig');
                    savefig(figName)
                else
                    warning(['Not all necessary tracking results present for ' filename ])
                end
            end
        end
    end
end