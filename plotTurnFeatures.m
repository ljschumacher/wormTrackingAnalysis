% plot features of red worms (from the second dataset) that may be useful
% for identifying omega turns

clear
close all

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

%% set parameters
phase = 'fullMovie'; % 'fullMovie' or 'stationary'. Script defines stationary phase as: starts at 10% into the movie, and stops at 60% into the movie (HA and N2) or at specified stopping frames (npr-1).
strains = {'npr1','N2'}; % {'npr1','N2'}
wormnums = {'40'};% {'40','HD'};

intensityThresholds = containers.Map({'40','HD','1W'},{60, 40, 100});
maxBlobSize = 2.5e5;
minSkelLength = 850;
maxSkelLength = 1500;
minNeighbrDist = 2000;
inClusterNeighbourNum = 3;
pixelsize = 100/19.5; % 100 microns are 19.5 pixels

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        [lastFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_r_list.xlsx'],1,'A1:B15','basic');
        numFiles = length(filenames);
        % create cell arrays to hold individual movie values to be pooled
        midbodyAngularSpeed_leaveCluster = cell(numFiles,1);
        pathCurvature_leaveCluster = cell(numFiles,1);
        midbodyAngularSpeed_loneWorms = cell(numFiles,1);
        pathCurvature_loneWorms = cell(numFiles,1);
        midbodyAngularSpeedFig = figure;
        pathCurvatureFig = figure;
        leaveClusterWormCount = 0;
        loneWormCount = 0;
       % go through individual movies
        for fileCtr = 1:numFiles
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton');
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            features = h5read(strrep(filename,'skeletons','features'),'/features_timeseries');
            if strcmp(phase, 'fullMovie')
                lastFrame = double(max(trajData.frame_number));
            elseif strcmp(phase,'stationary')
                lastFrame = lastFrames(fileCtr);
            end
            % filter red by blob size and intensity
            if contains(filename,'55')||contains(filename,'54')
                intensityThreshold = 80;
            else
                intensityThreshold = 40;
            end
            trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                intensityThreshold,maxBlobSize);
            % filter red by skeleton length
            trajData.filtered = trajData.filtered&logical(trajData.is_good_skel)...
                &filterSkelLength(skelData,pixelsize,minSkelLength,maxSkelLength);
            % restrict movies to stationary phase
            if strcmp(phase,'stationary')
                firstFrame = double(round(max(trajData.frame_number)/10)); % define starting frame as 10% into the movie
                phaseFrameLogInd = trajData.frame_number < lastFrame & trajData.frame_number > firstFrame;
                trajData.filtered(~phaseFrameLogInd) = false;
            end
            features.filtered = ismember(features.skeleton_id+1,find(trajData.filtered)); % use trajData.filtered to filter out unwa
            % restrict to worms that have just left a cluster
            min_neighbr_dist = h5read(filename,'/min_neighbr_dist');
            num_close_neighbrs = h5read(filename,'/num_close_neighbrs');
            neighbr_dist = h5read(filename,'/neighbr_distances');
            inClusterLogInd = num_close_neighbrs>=inClusterNeighbourNum;
            leaveClusterLogInd = vertcat(false,inClusterLogInd(1:end-1)&~inClusterLogInd(2:end)); % find worm-frames where inCluster changes from true to false
            leaveClusterFrameStart = find(leaveClusterLogInd); 
            leaveClusterFrameEnd = leaveClusterFrameStart+5*frameRate; %retain data for 5 seconds after a worm exits cluster
            leaveClusterFrameEnd = leaveClusterFrameEnd(leaveClusterFrameEnd<=numel(leaveClusterLogInd)); % exclude movie segments with frames beyond highest frame number
            leaveClusterFrameStart = leaveClusterFrameStart(1:numel(leaveClusterFrameEnd));
            for exitCtr = 1:numel(leaveClusterFrameStart)
                leaveClusterLogInd(leaveClusterFrameStart(exitCtr):leaveClusterFrameEnd(exitCtr))=true;
            end
            %the following line should do the same as the loop above, if dimensions of
            %terms added are compatible (eg row + column vector)
             %leaveClusterLogInd(unique(leaveClusterFrameStart + 0:5*frameRate)) = true;
            leaveClusterLogInd(inClusterLogInd)=false; % exclude when worms move back into a cluster 
            loneWormLogInd = min_neighbr_dist>=minNeighbrDist;
            leaveClusterLogInd(loneWormLogInd)=false; % exclude worms that have become lone worm
            leaveClusterLogInd = ismember(features.skeleton_id+1,find(trajData.filtered & leaveClusterLogInd)); % make clusterClusterLogInd the same size as features.filtered
            % restrict to lone worms
            loneWormLogInd = ismember(features.skeleton_id+1,find(trajData.filtered & loneWormLogInd));
            % write individual file results into cell array for pooling across movies later
            midbodyAngularSpeed_leaveCluster{fileCtr} = abs(features.midbody_motion_direction(features.filtered&leaveClusterLogInd));
            pathCurvature_leaveCluster{fileCtr} = abs(features.path_curvature(features.filtered&leaveClusterLogInd));
            midbodyAngularSpeed_loneWorms{fileCtr} = abs(features.midbody_motion_direction(features.filtered&loneWormLogInd));
            pathCurvature_loneWorms{fileCtr} = abs(features.path_curvature(features.filtered&loneWormLogInd));
            leaveClusterWormCount = leaveClusterWormCount + numel(unique(features.worm_index(leaveClusterLogInd)));
            loneWormCount = loneWormCount + numel(unique(features.worm_index(loneWormLogInd)));
        end
        % pool data from all files belonging to the same strain and worm density
        midbodyAngularSpeed_leaveCluster = vertcat(midbodyAngularSpeed_leaveCluster{:});
        pathCurvature_leaveCluster = vertcat(pathCurvature_leaveCluster{:});
        midbodyAngularSpeed_loneWorms = vertcat(midbodyAngularSpeed_loneWorms{:});
        pathCurvature_loneWorms = vertcat(pathCurvature_loneWorms{:});
        %% plot data
        % midbodyAngularSpeed figure
        set(0,'CurrentFigure',midbodyAngularSpeedFig)
        histogram(midbodyAngularSpeed_leaveCluster,'Normalization','pdf','DisplayStyle','stairs') 
        hold on
        histogram(midbodyAngularSpeed_loneWorms,'Normalization','pdf','DisplayStyle','stairs')
        leaveClusterLegend = ['leave cluster, n = ' num2str(leaveClusterWormCount)];
        loneWormLegend = ['lone worm, n = ' num2str(loneWormCount)];
        legend(leaveClusterLegend, loneWormLegend)
        title([strains{strainCtr} '\_' wormnums{numCtr} '\_midbodyAngularSpeed'],'FontWeight','normal')
        xlabel('midbody angular speed (degree/s)')
        ylabel('probability')
        xlim([0 8])
        ylim([0 2])
        set(midbodyAngularSpeedFig,'PaperUnits','centimeters')
        figurename = ['figures/turns/midbodyAngularSpeed_' strains{strainCtr} '_' wormnums{numCtr} '_' phase '_CL'];
        exportfig(midbodyAngularSpeedFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        
        % pathCurvature figure
        set(0,'CurrentFigure',pathCurvatureFig)
        histogram(pathCurvature_leaveCluster,'Normalization','pdf','DisplayStyle','stairs')
        hold on
        histogram(pathCurvature_loneWorms,'Normalization','pdf','DisplayStyle','stairs')
        legend(leaveClusterLegend, loneWormLegend)
        title([strains{strainCtr} '\_' wormnums{numCtr} '\_pathCurvature'],'FontWeight','normal')
        xlabel('path curvature (radians/microns)')
        ylabel('probability')
        xlim([0 0.5])
        ylim([0 60])
        set(pathCurvatureFig,'PaperUnits','centimeters')
        figurename = ['figures/turns/pathCurvature_' strains{strainCtr} '_' wormnums{numCtr} '_' phase '_CL'];
        exportfig(pathCurvatureFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        
        % T-test for equal means without assuming equal variance
        %[midbodyAngularSpeed_h, midbodyAngularSpeed_p] = ttest2(midbodyAngularSpeed_leaveCluster, midbodyAngularSpeed_loneWorms,'Vartype','unequal')
        %[pathCurvature_h, pathCurvature_p] = ttest2(pathCurvature_leaveCluster, pathCurvature_loneWorms,'Vartype','unequal')
    end
end