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
phase = 'fullMovie'; % 'fullMovie' or 'stationary'
strains = {'npr1'}; %{'npr1','N2'}
wormnums = {'40'};%{'40','HD'};

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
        % load data
        if strcmp(phase, 'stationary')
            [lastFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_r_list.xlsx'],1,'A1:B15','basic');
        elseif strcmp(phase,'fullMovie')
            filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_r_list.txt']);
        end
        numFiles = length(filenames);
        midbodyAngularSpeedFig = figure; hold on
        pathCurvatureFig = figure; hold on
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
                phaseFrameLogInd = trajData.frame_number < lastFrame;
                trajData.filtered(~phaseFrameLogInd) = false;
            end
            features.filtered = ismember(features.skeleton_id+1,find(trajData.filtered));
            % restrict to worms that have just left a cluster
            min_neighbr_dist = h5read(filename,'/min_neighbr_dist');
            num_close_neighbrs = h5read(filename,'/num_close_neighbrs');
            neighbr_dist = h5read(filename,'/neighbr_distances');
            loneWormLogInd = min_neighbr_dist>=minNeighbrDist;
            inClusterLogInd = num_close_neighbrs>=inClusterNeighbourNum;
            % find worm-frames where inCluster changes from true to false
            leaveClusterLogInd = vertcat(false,inClusterLogInd(1:end-1)&~inClusterLogInd(2:end));
            leaveClusterFrameStart = find(leaveClusterLogInd);
            leaveClusterFrameEnd = leaveClusterFrameStart+5*frameRate; %retain data for 5 seconds after a worm exits cluster
            for exitCtr = 1:numel(leaveClusterFrameStart)
                leaveClusterLogInd(leaveClusterFrameStart(exitCtr):leaveClusterFrameEnd(exitCtr))=true;
            end
%             the following line should do the same as the loop above, if dimensions of
%             terms added are compatible (eg row + column vector)
%             leaveClusterLogInd(unique(leaveClusterFrameStart + 0:5*frameRate)) = true;
            leaveClusterLogInd(inClusterLogInd)=false;
            % plot feature distributions
            set(0,'CurrentFigure',midbodyAngularSpeedFig)
            histogram(abs(features.midbody_motion_direction(features.filtered&leaveClusterLogInd)),'Normalization','pdf','DisplayStyle','stairs')
            set(0,'CurrentFigure',pathCurvatureFig)
            histogram(abs(features.path_curvature(features.filtered&leaveClusterLogInd)),'Normalization','pdf','DisplayStyle','stairs')
        end
%         % format and export figures
%         set(0,'CurrentFigure',midbodyAngularSpeedFig)
%         title([strains{strainCtr} '\_' wormnums{numCtr} '\_midbodyAngularSpeed'],'FontWeight','normal')
%         xlabel('midbody angular speed (degree/s)')
%         ylabel('probability')
%         xlim([0 8])
%         ylim([0 3])
%         set(midbodyAngularSpeedFig,'PaperUnits','centimeters')
%         figurename = ['figures/omegaTurns/midbodyAngularSpeed_' strains{strainCtr} '_' wormnums{numCtr} '_' phase];
%         exportfig(midbodyAngularSpeedFig,[figurename '.eps'],exportOptions)
%         system(['epstopdf ' figurename '.eps']);
%         system(['rm ' figurename '.eps']);
%         %
%         set(0,'CurrentFigure',pathCurvatureFig)
%         title([strains{strainCtr} '\_' wormnums{numCtr} '\_pathCurvature'],'FontWeight','normal')
%         xlabel('path curvature (radians/microns)')
%         ylabel('probability')
%         xlim([0 0.5])
%         ylim([0 140])
%         set(pathCurvatureFig,'PaperUnits','centimeters')
%         figurename = ['figures/omegaTurns/pathCurvature_' strains{strainCtr} '_' wormnums{numCtr} '_' phase];
%         exportfig(pathCurvatureFig,[figurename '.eps'],exportOptions)
%         system(['epstopdf ' figurename '.eps']);
%         system(['rm ' figurename '.eps']);
    end
end