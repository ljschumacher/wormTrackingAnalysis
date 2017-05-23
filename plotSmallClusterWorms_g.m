% plot small cluster worms with their neighbr
% and optionally, pharynx direction
% and/or their trajectories for the 100 framed preceding/following the randomly chosen frame


% remaining issue: when exporting as pdf via eps, some of the dots get
% lost, thus script is currently saving both .fig and .pdf files

close all
clear

% specify data sets
dataset = 2; % specify which dataset to run the script for. Enter either 1 or 2
strains = {'npr1','N2'};
wormnums = {'40','HD'};
numFramesSampled = 10; % how many frames to randomly sample per file

% set parameter for plotting trajectory/pharynx direction
plotPharynxDirection = true; % true or false
plotTraj = false; % true or false
trajFrameNumBefore = 100;
trajFrameNumAfter = 100;


% set parameters for filtering data
neighbrCutOff = 500; % distance in microns to consider a neighbr close
maxBlobSize = 1e4;
loneClusterRadius = 2000; % distance in microns to consider a cluster by itself
intensityThresholds = [60, 40];
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
maxBlobSize_r = 2.5e5;
minSkelLength_r = 850;
maxSkelLength_r = 1500;

% analysis
for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        if dataset == 1
            %% load data
            filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_list.txt']);
        elseif dataset == 2
            %% load data
            filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_g_list.txt']);
            filenames_r = importdata(['datalists/' strains{strainCtr} '_' wormnum '_r_list.txt']);
        end
        numFiles = length(filenames);
        for fileCtr=1:numFiles
            if dataset ==1
                filename = filenames{fileCtr};
                trajData = h5read(filename,'/trajectories_data');
                skelData = h5read(filename,'/skeleton');
                blobFeats = h5read(filename,'/blob_features');
                numCloseNeighbr = h5read(filename,'/num_close_neighbrs');
                neighbrDist = h5read(filename,'/neighbr_distances');
            elseif dataset ==2
                filename = filenames{fileCtr};
                filename_r = filenames_r{fileCtr};
                if exist(filename,'file')&&exist(filename_r,'file')
                    trajData = h5read(filename,'/trajectories_data');
                    skelData = h5read(filename,'/skeleton');
                    blobFeats = h5read(filename,'/blob_features');
                    numCloseNeighbr = h5read(filename,'/num_close_neighbrs');
                    neighbrDist = h5read(filename,'/neighbr_distances');
                    trajData_r = h5read(filename_r,'/trajectories_data');
                    blobFeats_r = h5read(filename_r,'/blob_features');
                    skelData_r = h5read(filename_r,'/skeleton');
                end
            end
            %% filter data
            % filter green by blob size and intensity
            trajData.filtered = (blobFeats.area*pixelsize^2<=maxBlobSize)&...
                (blobFeats.intensity_mean>=intensityThresholds(numCtr));
            % filter green by small cluster status
            trajData.clusterfiltered = trajData.filtered&...
                ((numCloseNeighbr== 2 & neighbrDist(:,3)>=loneClusterRadius)...
                |(numCloseNeighbr== 3 & neighbrDist(:,4)>=(loneClusterRadius))...
                |(numCloseNeighbr== 4 & neighbrDist(:,5)>=(loneClusterRadius)));
            if dataset ==2
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
            end
            % plot sample data
            if length(unique(trajData.frame_number(trajData.clusterfiltered)))<numFramesSampled
                warning(['Not enough frames to plot for ' filename ])
            else
                smallClusterTrajFig = figure;
                framesAnalyzed = randsample(unique(trajData.frame_number(trajData.clusterfiltered)),numFramesSampled);
                for frameCtr = 1:numFramesSampled
                    frame=framesAnalyzed(frameCtr);
                    subplot(2,5,frameCtr)
                    % plot central green worm (in blue)
                    frameIdcs_worm = find(trajData.frame_number==frame&trajData.clusterfiltered);
                    if nnz(frameIdcs_worm)>1 % plot only one green worm if multiple present
                        frameIdcs_worm = randsample(frameIdcs_worm,1);
                    end
                    worm_xcoord = trajData.coord_x(frameIdcs_worm);
                    worm_ycoord = trajData.coord_y(frameIdcs_worm);
                    if plotPharynxDirection == false
                        plot(worm_xcoord,worm_ycoord,'ko','MarkerSize',5,'MarkerFaceColor','b')
                        hold on
                    elseif plotPharynxDirection == true
                        if isnan(skelData(1,1,frameIdcs_worm)) == true || isnan(skelData(2,1,frameIdcs_worm)) == true
                            warning('worm not plotted!')
                        else
                        plot(skelData(1,1,frameIdcs_worm),skelData(2,1,frameIdcs_worm),'ko','MarkerSize',3,'MarkerFaceColor','b') %plot head
                        hold on
                        plot(skelData(1,:,frameIdcs_worm),skelData(2,:,frameIdcs_worm),'LineWidth',2,'Color','b')%plot pharynx
                        end
                    end
                    % plot central worm trajectory
                    if plotTraj == true
                        centTrajIdcsBefore = (frameIdcs_worm-trajFrameNumBefore):frameIdcs_worm;
                        centTrajIdcsBefore = centTrajIdcsBefore(centTrajIdcsBefore>0);
                        centTrajIdcsAfter = frameIdcs_worm:(frameIdcs_worm+trajFrameNumAfter);
                        sameCentWormBefore = trajData.worm_index_joined(centTrajIdcsBefore);
                        centTrajIdcsBefore(sameCentWormBefore ~=trajData.worm_index_joined(frameIdcs_worm))=[];
                        sameWormAfter = trajData.worm_index_joined(centTrajIdcsAfter);
                        centTrajIdcsAfter(sameWormAfter ~=trajData.worm_index_joined(frameIdcs_worm))=[];
                        % only keep trajectory indices for the same worm
                        plot(trajData.coord_x(centTrajIdcsBefore),trajData.coord_y(centTrajIdcsBefore),...
                            'Color',[0.4 0.4 1],'LineWidth',1)
                        plot(trajData.coord_x(centTrajIdcsAfter),trajData.coord_y(centTrajIdcsAfter),...
                            'Color','black','LineWidth',1)
                    end
                    % plot other green worms
                    frameLogIdcs_pharynx = trajData.frame_number==frame&trajData.filtered;
                    frameLogIdcs_pharynx(frameIdcs_worm) = false; % exclude the central worm that's just been plotted
                    if plotPharynxDirection == false
                        plot(trajData.coord_x(frameLogIdcs_pharynx),trajData.coord_y(frameLogIdcs_pharynx),...
                            'ro','MarkerSize',5,'MarkerFaceColor',[0 0.7 0.3])
                    elseif plotPharynxDirection == true
                        head_x = squeeze(skelData(1,1,frameLogIdcs_pharynx));
                        head_y = squeeze(skelData(2,1,frameLogIdcs_pharynx));
                        plot(head_x,head_y,'ko','MarkerSize',3,'MarkerFaceColor',[0 0.7 0.3]) %plot head
                        pharynxSkel_x = squeeze(skelData(1,:,frameLogIdcs_pharynx));
                        pharynxSkel_y = squeeze(skelData(2,:,frameLogIdcs_pharynx));
                        plot(pharynxSkel_x,pharynxSkel_y,'LineWidth',2,'Color',[0 0.7 0.3])%plot pharynx
                    end
                    % plot other green worm trajectories
                    if plotTraj == true
                        greenTrajIdcsList = find(frameLogIdcs_pharynx);
                        numGreenTrajIdcs = length(greenTrajIdcsList);
                        for greenTrajCtr = 1:numGreenTrajIdcs % loop through each green worm
                            greenTrajInd = greenTrajIdcsList(greenTrajCtr);
                            greenTrajIdcsBefore = (greenTrajInd-trajFrameNumBefore):greenTrajInd;
                            greenTrajIdcsBefore = greenTrajIdcsBefore(greenTrajIdcsBefore>0);
                            sameGreenWormBefore = trajData.worm_index_joined(greenTrajIdcsBefore);
                            greenTrajIdcsBefore(sameGreenWormBefore ~=trajData.worm_index_joined(greenTrajInd))=[];
                            greenTrajIdcsAfter = greenTrajInd:(greenTrajInd+trajFrameNumAfter);
                            sameGreenWormAfter = trajData.worm_index_joined(greenTrajIdcsAfter);
                            greenTrajIdcsAfter(sameGreenWormAfter ~=trajData.worm_index_joined(greenTrajInd))=[];
                            % only keep trajectory indices for the same worm
                            plot(trajData.coord_x(greenTrajIdcsBefore),trajData.coord_y(greenTrajIdcsBefore),...
                                'Color','green','LineWidth',1)
                            plot(trajData.coord_x(greenTrajIdcsAfter),trajData.coord_y(greenTrajIdcsAfter),...
                                'Color',[0.2 0.6 0.2],'LineWidth',1)
                        end
                    end
                    if dataset ==2
                        % plot red worms
                        frameLogIdcs_red = trajData_r.frame_number==frame&trajData_r.filtered;
                        plot(trajData_r.coord_x(frameLogIdcs_red),trajData_r.coord_y(frameLogIdcs_red),...
                            'ko','MarkerSize',4,'MarkerFaceColor','m')
                        % plot red worm trajectories
                        if plotTraj == true
                            redTrajIdcsList = find(frameLogIdcs_red);
                            numRedTrajIdcs = length(redTrajIdcsList);
                            for redTrajCtr = 1:numRedTrajIdcs % loop through each red worm
                                redTrajInd = redTrajIdcsList(redTrajCtr);
                                redTrajIdcsBefore = (redTrajInd-trajFrameNumBefore):redTrajInd;
                                redTrajIdcsBefore = redTrajIdcsBefore(redTrajIdcsBefore>0);
                                sameRedWormBefore = trajData.worm_index_joined(redTrajIdcsBefore);
                                redTrajIdcsBefore(sameRedWormBefore ~=trajData.worm_index_joined(redTrajInd))=[];
                                redTrajIdcsAfter = redTrajInd:(redTrajInd+trajFrameNumAfter);
                                sameRedWormAfter = trajData.worm_index_joined(redTrajIdcsAfter);
                                redTrajIdcsAfter(sameRedWormAfter ~=trajData.worm_index_joined(redTrajInd))=[];
                                % only keep trajectory indices for the same worm
                                plot(trajData_r.coord_x(redTrajIdcsBefore),trajData_r.coord_y(redTrajIdcsBefore),...
                                    'Color',[1 0 0.2],'LineWidth',1)
                                plot(trajData_r.coord_x(redTrajIdcsAfter),trajData_r.coord_y(redTrajIdcsAfter),...
                                    'Color',[0.6 0 0.2],'LineWidth',1)
                            end
                        end
                    end
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
                set(smallClusterTrajFig,'Name',[strain ' ' wormnum ' ' figName])
                if plotTraj == true
                    if dataset ==1
                        epsFileName = ['figures/smallClusterTraj/green1/pdf/sampleSmallClusterTraj_' strain '_' wormnum '_' figName '.eps'];
                        figFileName = ['figures/smallClusterTraj/green1/fig/sampleSmallClusterTraj_' strain '_' wormnum '_' figName '.fig'];
                    elseif dataset ==2
                        epsFileName = ['figures/smallClusterTraj/green2/pdf/sampleSmallClusterTraj_' strain '_' wormnum '_' figName '.eps'];
                        figFileName = ['figures/smallClusterTraj/green2/fig/sampleSmallClusterTraj_' strain '_' wormnum '_' figName '.fig'];
                    end
                elseif plotTraj == false
                    if dataset ==1
                        epsFileName = ['figures/smallCluster/green1/pdf/sampleSmallCluster_' strain '_' wormnum '_' figName '.eps'];
                        figFileName = ['figures/smallCluster/green1/fig/sampleSmallCluster_' strain '_' wormnum '_' figName '.fig'];
                    elseif dataset ==2
                        epsFileName = ['figures/smallCluster/green2/pdf/sampleSmallCluster_' strain '_' wormnum '_' figName '.eps'];
                        figFileName = ['figures/smallCluster/green2/fig/sampleSmallCluster_' strain '_' wormnum '_' figName '.fig'];
                    end
                end
                %savefig(figFileName)
                %exportfig(smallClusterTrajFig,epsFileName,'Color','rgb','Width',210,'Height',297)
                %system(['epstopdf ' epsFileName]);
                %system(['rm ' epsFileName]);
                %close all
            end
        end
    end
end