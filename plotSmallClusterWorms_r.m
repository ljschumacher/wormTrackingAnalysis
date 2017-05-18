% plot small cluster worms with their neighbrs
% and optionally, their trajectories for the 100 framed preceding/following the randomly chosen frame

% remaining issue: when exporting as pdf via eps, some of the dots get
% lost, thus script is currently saving both .fig and .pdf files

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

% set parameter for plotting trajectory
plotTraj = true; % true or false
trajFrameNumBefore = 100;
trajFrameNumAfter = 100;

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        %% load data
        filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_r_list.txt']);
        numFiles = length(filenames);
        filenames_g = importdata(['datalists/' strains{strainCtr} '_' wormnum '_g_list.txt']);
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
                % plot sample data
                if length(unique(trajData.frame_number(trajData.filtered)))<numFramesSampled
                    warning(['Not enough frames to plot for ' filename ])
                else
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
                        plot(worm_xcoords,worm_ycoords,'Color',[0.8 0 0.2],'LineWidth',3)
                        hold on
                        % plot red worm trajectory
                        if plotTraj == true;
                            redTrajIdcsBefore = (frameIdcs_worm-trajFrameNumBefore):frameIdcs_worm;
                            redTrajIdcsBefore = redTrajIdcsBefore(redTrajIdcsBefore>0); %remove negative values
                            redTrajIdcsAfter = frameIdcs_worm:(frameIdcs_worm+trajFrameNumAfter);
                            sameRedWormBefore = trajData.worm_index_joined(redTrajIdcsBefore);
                            redTrajIdcsBefore(sameRedWormBefore ~=trajData.worm_index_joined(frameIdcs_worm))=[];
                            sameWormAfter = trajData.worm_index_joined(redTrajIdcsAfter);
                            redTrajIdcsAfter(sameWormAfter ~=trajData.worm_index_joined(frameIdcs_worm))=[];
                            % only keep trajectory indices for the same worm
                            plot(trajData.coord_x(redTrajIdcsBefore),trajData.coord_y(redTrajIdcsBefore),...
                                'Color',[1 0 0.2],'LineWidth',1)
                            plot(trajData.coord_x(redTrajIdcsAfter),trajData.coord_y(redTrajIdcsAfter),...
                                'Color',[0.6 0 0.2],'LineWidth',1)
                            plot(trajData.coord_x(frameIdcs_worm),trajData.coord_y(frameIdcs_worm),...
                                'o','MarkerSize',8,'MarkerFaceColor',[0.8 0 0.2])
                        end
                        % plot green worms
                        frameLogIdcs_pharynx = trajData_g.frame_number==frame&trajData_g.filtered;
                        plot(trajData_g.coord_x(frameLogIdcs_pharynx),trajData_g.coord_y(frameLogIdcs_pharynx),...
                            'ko','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','k')
                        % plot green worm trajectories
                        if plotTraj == true;
                            greenTrajIdcsList = find(frameLogIdcs_pharynx);
                            numGreenTrajIdcs = length(greenTrajIdcsList);
                            for greenTrajCtr = 1:numGreenTrajIdcs % loop through each green worm
                                greenTrajInd = greenTrajIdcsList(greenTrajCtr);
                                greenTrajIdcsBefore = (greenTrajInd-trajFrameNumBefore):greenTrajInd;
                                greenTrajIdcsBefore = greenTrajIdcsBefore(greenTrajIdcsBefore>0); %remove negative values
                                sameGreenWormBefore = trajData_g.worm_index_joined(greenTrajIdcsBefore);
                                greenTrajIdcsBefore(sameGreenWormBefore ~=trajData_g.worm_index_joined(greenTrajInd))=[];
                                greenTrajIdcsAfter = greenTrajInd:(greenTrajInd+trajFrameNumAfter);
                                sameGreenWormAfter = trajData_g.worm_index_joined(greenTrajIdcsAfter);
                                greenTrajIdcsAfter(sameGreenWormAfter ~=trajData_g.worm_index_joined(greenTrajInd))=[];
                                % only keep trajectory indices for the same worm
                                plot(trajData_g.coord_x(greenTrajIdcsBefore),trajData_g.coord_y(greenTrajIdcsBefore),...
                                    'Color',[0.2 0 0.8],'LineWidth',1)
                                plot(trajData_g.coord_x(greenTrajIdcsAfter),trajData_g.coord_y(greenTrajIdcsAfter),...
                                    'Color',[0.2 0 0.2],'LineWidth',1)
                            end
                        end
                        axis equal
                        % plot circle of radius neighbrCutOff around each worm
                        viscircles([trajData.coord_x(frameIdcs_worm) trajData.coord_y(frameIdcs_worm)],...
                            neighbrCutOff/pixelsize,'LineStyle','--','Color',0.5*[1 1 1],'EnhanceVisibility',false);
                        % plot circle of radius loneClusterRadius around each worm
                        viscircles([trajData.coord_x(frameIdcs_worm) trajData.coord_y(frameIdcs_worm)],...
                            loneClusterRadius/pixelsize,'LineStyle','--','Color',0.5*[1 1 1],'EnhanceVisibility',false);
                        % %%% plot format
                        ax = gca;
                        xlim([-2000 2000]/pixelsize + trajData.coord_x(frameIdcs_worm))
                        ylim([-2000 2000]/pixelsize + trajData.coord_y(frameIdcs_worm))
                        set(ax,'visible','off')
                        ax.Position = ax.Position.*[1 1 1.2 1.2]; % reduce whitespace btw subplots
                    end
                    %% export figure
                    figName = strrep(strrep(filename(end-32:end-17),'_',''),'/','');
                    set(smallClusterWormsFig,'Name',[strain ' ' wormnum ' ' figName])
                    if plotTraj == true
                        epsFileName = ['figures/smallClusterTraj/red2/pdf/sampleSmallClusterTraj_' strain '_' wormnum '_' figName '.eps'];
                        figFileName = ['figures/smallClusterTraj/red2/fig/sampleSmallClusterTraj_' strain '_' wormnum '_' figName '.fig'];
                    elseif plotTraj == false
                        epsFileName = ['figures/smallCluster/red2/pdf/sampleSmallCluster_' strain '_' wormnum '_' figName '.eps'];
                        figFileName = ['figures/smallCluster/red2/fig/sampleSmallCluster_' strain '_' wormnum '_' figName '.fig'];
                    end
                    savefig(figFileName)
                    exportfig(smallClusterWormsFig,epsFileName,'Color','rgb','Width',210,'Height',297)
                    system(['epstopdf ' epsFileName]);
                    system(['rm ' epsFileName]);
                end
            end
        end
    end
end