% plot small cluster worms with their neighbrs
close all
clear

% specify data sets
dataset = 1 % specify which dataset to run the script for. Enter either 1 or 2
strains = {'npr1','N2'};
wormnums = {'40','HD'};
numFramesSampled = 30; % how many frames to randomly sample per file

% set parameters for filtering data
neighbrCutOff = 500; % distance in microns to consider a neighbr close
maxBlobSize = 1e4;
loneClusterRadius = 2000; % distance in microns to consider a cluster by itself
intensityThresholds = [60, 40];
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
maxBlobSize_r = 2.5e5;
minSkelLength_r = 850;
maxSkelLength_r = 1500;

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        if dataset == 1
            %% load data
            filenames = importdata([strains{strainCtr} '_' wormnum '_list.txt']);
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
                trajData.clusterfiltered = trajData.filtered&...
                    ((numCloseNeighbr== 2 & neighbrDist(:,3)>=loneClusterRadius)...
                    |(numCloseNeighbr== 3 & neighbrDist(:,4)>=(loneClusterRadius))...
                    |(numCloseNeighbr== 4 & neighbrDist(:,5)>=(loneClusterRadius)));
                % plot sample data
                if length(unique(trajData.frame_number(trajData.clusterfiltered)))<numFramesSampled
                    warning(['Not enough frames to plot for ' filename ])
                else
                    smallClusterWormsFig = figure;
                    framesAnalyzed = randsample(unique(trajData.frame_number(trajData.clusterfiltered)),numFramesSampled);
                    for frameCtr = 1:numFramesSampled
                        frame=framesAnalyzed(frameCtr);
                        subplot(6,5,frameCtr)
                        % plot central green worm (in blue)
                        frameIdcs_worm = find(trajData.frame_number==frame&trajData.clusterfiltered);
                        if nnz(frameIdcs_worm)>1 % plot only one green worm if multiple present
                            frameIdcs_worm = randsample(frameIdcs_worm,1);
                        end
                        worm_xcoord = trajData.coord_x(frameIdcs_worm);
                        worm_ycoord = trajData.coord_y(frameIdcs_worm);
                        plot(worm_xcoord,worm_ycoord,'ko','MarkerSize',4,'MarkerFaceColor','b')
                        hold on
                        % plot other green worms
                        frameLogIdcs_pharynx = trajData.frame_number==frame&trajData.filtered;
                        frameLogIdcs_pharynx(frameIdcs_worm) = false; % exclude the central worm that's just been plotted
                        plot(trajData.coord_x(frameLogIdcs_pharynx),trajData.coord_y(frameLogIdcs_pharynx),...
                            'ko','MarkerSize',4,'MarkerFaceColor','g')
                        axis equal
                        % plot circle of radius neighbrCutOff around each worm
                        viscircles([worm_xcoord worm_ycoord],...
                            neighbrCutOff/pixelsize,'LineStyle','--','Color',0.5*[1 1 1],'EnhanceVisibility',false);
                        % plot circle of radius loneClusterRadius around each worm
                        viscircles([worm_xcoord worm_ycoord],...
                            loneClusterRadius/pixelsize,'LineStyle','--','Color',0.5*[1 1 1],'EnhanceVisibility',false);
                        % %%% plot format
                        ax = gca;
                        xlim([-2000 2000]/pixelsize + worm_xcoord)
                        ylim([-2000 2000]/pixelsize + worm_ycoord)
                        set(ax,'visible','off')
                        ax.Position = ax.Position.*[1 1 1.2 1.2]; % reduce whitespace btw subplots
                    end
                    %% export figure
                    figName = strrep(strrep(filename(end-32:end-17),'_',''),'/','');
                    set(smallClusterWormsFig,'Name',[strain ' ' wormnum ' ' figName])
                    epsFileName = ['figures/diagnostics/smallCluster/green1/pdf/sampleSmallClusterWorms_' strain '_' wormnum '_' figName '.eps'];
                    figFileName = ['figures/diagnostics/smallCluster/green1/fig/sampleSmallClusterWorms_' strain '_' wormnum '_' figName '.fig'];
                    savefig(figFileName)
                    exportfig(smallClusterWormsFig,epsFileName,'Color','rgb','Width',210,'Height',297)
                    system(['epstopdf ' epsFileName]);
                    system(['rm ' epsFileName]);
                end
            end
        elseif dataset == 2
            %% load data
            filenames = importdata([strains{strainCtr} '_' wormnum '_g_list.txt']);
            filenames_r = importdata([strains{strainCtr} '_' wormnum '_r_list.txt']);
            numFiles = length(filenames);
            for fileCtr=1:numFiles
                filename = filenames{fileCtr};
                filename_r = filenames_r{fileCtr};
                if exist(filename,'file')&&exist(filename_r,'file')
                    trajData = h5read(filename,'/trajectories_data');
                    blobFeats = h5read(filename,'/blob_features');
                    numCloseNeighbr = h5read(filename,'/num_close_neighbrs');
                    neighbrDist = h5read(filename,'/neighbr_distances');
                    trajData_r = h5read(filename_r,'/trajectories_data');
                    blobFeats_r = h5read(filename_r,'/blob_features');
                    skelData_r = h5read(filename_r,'/skeleton');
                    %% filter data
                    % filter green by blob size and intensity
                    trajData.filtered = (blobFeats.area*pixelsize^2<=maxBlobSize)&...
                        (blobFeats.intensity_mean>=intensityThresholds(numCtr));
                    % filter green by small cluster status
                    trajData.clusterfiltered = trajData.filtered&...
                        ((numCloseNeighbr== 2 & neighbrDist(:,3)>=loneClusterRadius)...
                        |(numCloseNeighbr== 3 & neighbrDist(:,4)>=(loneClusterRadius))...
                        |(numCloseNeighbr== 4 & neighbrDist(:,5)>=(loneClusterRadius)));
                    % filter red by blob size and intensity
                    if contains(filename,'55')||contains(filename,'54')
                        intensityThreshold_r = 80;
                    else
                        intensityThreshold_r = 40;
                    end
                    trajData_r.filtered = filterIntensityAndSize(blobFeats_r,pixelsize,...
                        intensityThreshold_r,maxBlobSize_r);
                    % filter red by skeleton length
                    trajData_r.filtered = trajData_r.filtered&logical(trajData_r.is_good_skel)...
                        &filterSkelLength(skelData_r,pixelsize,minSkelLength_r,maxSkelLength_r);
                    % plot sample data
                    if length(unique(trajData.frame_number(trajData.clusterfiltered)))<numFramesSampled
                        warning(['Not enough frames to plot for ' filename ])
                    else
                        smallClusterWormsFig = figure;
                        framesAnalyzed = randsample(unique(trajData.frame_number(trajData.clusterfiltered)),numFramesSampled);
                        for frameCtr = 1:numFramesSampled
                            frame=framesAnalyzed(frameCtr);
                            subplot(6,5,frameCtr)
                            % plot central green worm (in blue)
                            frameIdcs_worm = find(trajData.frame_number==frame&trajData.clusterfiltered);
                            if nnz(frameIdcs_worm)>1 % plot only one green worm if multiple present
                                frameIdcs_worm = randsample(frameIdcs_worm,1);
                            end
                            worm_xcoord = trajData.coord_x(frameIdcs_worm);
                            worm_ycoord = trajData.coord_y(frameIdcs_worm);
                            plot(worm_xcoord,worm_ycoord,'ko','MarkerSize',4,'MarkerFaceColor','b')
                            hold on
                            % plot other green worms
                            frameLogIdcs_pharynx = trajData.frame_number==frame&trajData.filtered;
                            frameLogIdcs_pharynx(frameIdcs_worm) = false; % exclude the central worm that's just been plotted
                            plot(trajData.coord_x(frameLogIdcs_pharynx),trajData.coord_y(frameLogIdcs_pharynx),...
                                'ko','MarkerSize',4,'MarkerFaceColor','g')
                            % plot red worms
                            frameLogIdcs_red = trajData_r.frame_number==frame&trajData_r.filtered;
                            plot(trajData_r.coord_x(frameLogIdcs_red),trajData_r.coord_y(frameLogIdcs_red),...
                                'ko','MarkerSize',4,'MarkerFaceColor','m')
                            axis equal
                            % plot circle of radius neighbrCutOff around each worm
                            viscircles([worm_xcoord worm_ycoord],...
                                neighbrCutOff/pixelsize,'LineStyle','--','Color',0.5*[1 1 1],'EnhanceVisibility',false);
                            % plot circle of radius loneClusterRadius around each worm
                            viscircles([worm_xcoord worm_ycoord],...
                                loneClusterRadius/pixelsize,'LineStyle','--','Color',0.5*[1 1 1],'EnhanceVisibility',false);
                            % %%% plot format
                            ax = gca;
                            xlim([-2000 2000]/pixelsize + worm_xcoord)
                            ylim([-2000 2000]/pixelsize + worm_ycoord)
                            set(ax,'visible','off')
                            ax.Position = ax.Position.*[1 1 1.2 1.2]; % reduce whitespace btw subplots
                        end
                        %% export figure
                        figName = strrep(strrep(filename(end-32:end-17),'_',''),'/','');
                        set(smallClusterWormsFig,'Name',[strain ' ' wormnum ' ' figName])
                        epsFileName = ['figures/diagnostics/smallCluster/green2/pdf/sampleSmallClusterWorms_' strain '_' wormnum '_' figName '.eps'];
                        figFileName = ['figures/diagnostics/smallCluster/green2/fig/sampleSmallClusterWorms_' strain '_' wormnum '_' figName '.fig'];
                        savefig(figFileName)
                        exportfig(smallClusterWormsFig,epsFileName,'Color','rgb','Width',210,'Height',297)
                        system(['epstopdf ' epsFileName]);
                        system(['rm ' epsFileName]);
                    end
                end
            end
        end
    end
end