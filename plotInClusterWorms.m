% plot in cluster worms
close all
clear

% specify data sets
strains = {'npr1','N2'};
wormnums = {'40','HD'};
numFramesSampled = 30; % how many frames to randomly sample per file

% set parameters for filtering data
neighbrCutOff = 500; % distance in microns to consider a neighbr close
minNumNeighbrs = 3;
maxBlobSize = 2.5e5;
maxBlobSize_g = 1e4;
minSkelLength = 850;
maxSkelLength = 1500;
intensityThresholds_g = containers.Map({'40','HD','1W'},{60, 40, 100});
pixelsize = 100/19.5; % 100 microns are 19.5 pixels

for wormnum = wormnums
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};     
        %% load data
        filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum{1} '_r_list.txt']);
        numFiles = length(filenames);
        filenames_g = importdata(['datalists/' strains{strainCtr} '_' wormnum{1} '_g_list.txt']);
        for fileCtr=1:numFiles
            filename = filenames{fileCtr};
            filename_g = filenames_g{fileCtr};
            inClusterWormsFig = figure;
            if exist(filename,'file')&&exist(filename_g,'file')
                trajData = h5read(filename,'/trajectories_data');
                blobFeats = h5read(filename,'/blob_features');
                skelData = h5read(filename,'/skeleton');
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
                % filter green channel by blob size and intensity
                trajData_g.filtered = (blobFeats_g.area*pixelsize^2<=maxBlobSize_g)&...
                    (blobFeats_g.intensity_mean>=intensityThresholds_g(wormnum{1}));
                % filter for in-cluster
                num_close_neighbrs_rg = h5read(filename,'/num_close_neighbrs_rg')';
                trajData.filtered = trajData.filtered&num_close_neighbrs_rg>=minNumNeighbrs;
                % plot sample data
                framesAnalyzed = randsample(unique(trajData.frame_number(trajData.filtered)),numFramesSampled);
                for frameCtr = 1:numFramesSampled
                    frame=framesAnalyzed(frameCtr);
                    subplot(5,6,frameCtr)
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
                    %%% make separate plots for worms in clusters / lone
                    ax = gca;
                    xlim([-500 500]/pixelsize + trajData.coord_x(frameIdcs_worm))
                    ylim([-500 500]/pixelsize + trajData.coord_y(frameIdcs_worm))
                    set(ax,'visible','off')
                    ax.Position = ax.Position.*[1 1 1.2 1.2]; % reduce whitespace btw subplots
                end
                %% export figure
                figName = strrep(strrep(filename(end-32:end-17),'_',''),'/','');
                set(inClusterWormsFig,'Name',[strain ' ' wormnum{1} ' ' figName])
                figFileName = ['figures/diagnostics/sampleInClusterWorms_' strain '_' wormnum{1} '_' figName '.eps'];
                exportfig(inClusterWormsFig,figFileName,'Color','rgb')
                system(['epstopdf ' figFileName]);
                system(['rm ' figFileName]);
            else
                warning(['Not all necessary tracking results present for ' filename ])
            end
        end
    end
end
tilefigs()