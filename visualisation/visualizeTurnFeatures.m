% plot features of red worms (from the second dataset) that may be useful
% for identifying omega turns

% problem is that omega turn events from features file aren't time series
% thus cannot correlate with worm category restrictions. i.e. worm with
% index 1's overall omega turn frequency over its entire trajectory can be
% plotted, but information of worm index 1 during leave cluster or lone
% worm periods cannot be isolated, thus rendering omega turn feature less
% useful for our purpose. 

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
phase = 'fullMovie'; % 'fullMovie', 'joining', or 'sweeping'.
strains = {'npr1','N2'}; % {'npr1','N2'}
wormnums = {'40'};% {'40','HD'};
postExitDuration = 5; % set the duration (in seconds) after a worm exits a cluster to be included in the analysis

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
        [phaseFrames,filenames,~] = xlsread(['../datalists/' strains{strainCtr} '_' wormnum '_r_list.xlsx'],1,'A1:E15','basic');
        numFiles = length(filenames);
        % create cell arrays to hold individual movie values to be pooled
        maxPathNum = 200; % assume maximum number of paths per file is as set. Script checks for this later and gives warning if this is not enough
        headAngularSpeed_leaveCluster = cell(numFiles,1);
        headAngularSpeed_loneWorms = cell(numFiles,1);
        omegaTurnFreq_leaveCluster = cell(numFiles,1);
        omegaTurnFreq_loneWorm = cell(numFiles,1);
        headAngularSpeedFig = figure;
        omegaTurnFreqFig = figure;
        leaveClusterWormCount = 0;
        loneWormCount = 0;
        %% go through individual movies
        for fileCtr = 1:numFiles
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton');
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            features = h5read(strrep(filename,'skeletons','features'),'/features_timeseries');
            %% filter worms by various criteria
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
            % apply phase restriction
            [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr);
            phaseFrameLogInd = trajData.frame_number <= lastFrame & trajData.frame_number >= firstFrame;
            trajData.filtered(~phaseFrameLogInd) = false;
            features.filtered = ismember(features.skeleton_id+1,find(trajData.filtered)); % use trajData.filtered to filter features file
            % find worms that have just left a cluster
            [leaveClusterLogInd, loneWormLogInd, ~,~] = findWormCategory(filename,inClusterNeighbourNum,minNeighbrDist,postExitDuration);
            % turn the LogInd to the same size as features.filtered
            leaveClusterLogInd = ismember(features.skeleton_id+1,find(trajData.filtered & leaveClusterLogInd)); 
            loneWormLogInd = ismember(features.skeleton_id+1,find(trajData.filtered & loneWormLogInd));
            %% calculate or extract desired feature values
            % write additional turn features into cell array for pooling across movies later
            headAngularSpeed_leaveCluster{fileCtr} = abs(features.head_motion_direction(features.filtered&leaveClusterLogInd));
            headAngularSpeed_loneWorms{fileCtr} = abs(features.head_motion_direction(features.filtered&loneWormLogInd));
            uniqueLeaveClusterWorm = unique(features.worm_index(features.filtered&leaveClusterLogInd));
            uniqueLoneWorm = unique(features.worm_index(features.filtered&loneWormLogInd));
            omegaTurnFreq_leaveCluster_thisFile = zeros(1,numel(uniqueLeaveClusterWorm));
            for wormCtr = 1:numel(uniqueLeaveClusterWorm)
                omegaTurnFreq_leaveCluster_thisFile(wormCtr)= h5read(strrep(filename,'skeletons','features'),['/features_events/worm_' num2str(uniqueLeaveClusterWorm(wormCtr)) '/omega_turns_frequency']);
            end
            omegaTurnFreq_leaveCluster{fileCtr} = omegaTurnFreq_leaveCluster_thisFile;
            omegaTurnFreq_loneWorm_thisFile = zeros(1,numel(uniqueLoneWorm));
            for wormCtr = 1:numel(uniqueLoneWorm)
                omegaTurnFreq_loneWorm_thisFile(wormCtr)= h5read(strrep(filename,'skeletons','features'),['/features_events/worm_' num2str(uniqueLoneWorm(wormCtr)) '/omega_turns_frequency']);
            end
            omegaTurnFreq_loneWorm{fileCtr} = omegaTurnFreq_loneWorm_thisFile;
            % add number of worms for total n number across movies
            leaveClusterWormCount = leaveClusterWormCount + numel(unique(features.worm_index(leaveClusterLogInd)));
            loneWormCount = loneWormCount + numel(unique(features.worm_index(loneWormLogInd)));
        end
        %% pool data from all files belonging to the same strain and worm density
        headAngularSpeed_leaveCluster = vertcat(headAngularSpeed_leaveCluster{:});
        headAngularSpeed_loneWorms = vertcat(headAngularSpeed_loneWorms{:});
        omegaTurnFreq_leaveCluster = horzcat(omegaTurnFreq_leaveCluster{:})';
        omegaTurnFreq_loneWorm = horzcat(omegaTurnFreq_loneWorm{:})';
        %% plot data, format, and export
        %headAngularSpeed figure
        set(0,'CurrentFigure',headAngularSpeedFig)
        histogram(headAngularSpeed_leaveCluster,'Normalization','pdf','DisplayStyle','stairs')
        hold on
        histogram(headAngularSpeed_loneWorms,'Normalization','pdf','DisplayStyle','stairs')
        leaveClusterLegend = ['leave cluster, n = ' num2str(leaveClusterWormCount)];
        loneWormLegend = ['lone worm, n = ' num2str(loneWormCount)];
        legend(leaveClusterLegend, loneWormLegend)
        title([strains{strainCtr} '\_' wormnums{numCtr} '\_headAngularSpeed'],'FontWeight','normal')
        xlabel('head angular speed (degree/s)')
        ylabel('probability')
        xlim([0 8])
        ylim([0 2])
        set(headAngularSpeedFig,'PaperUnits','centimeters')
        figurename = ['../figures/turns/headAngularSpeed_' strains{strainCtr} '_' wormnums{numCtr} '_' phase '_CL'];
        %exportfig(headAngularSpeedFig,[figurename '.eps'],exportOptions)
        %system(['epstopdf ' figurename '.eps']);
        %system(['rm ' figurename '.eps']);
        
        %omegaTurns figure
        set(0,'CurrentFigure',omegaTurnFreqFig)
        histogram(omegaTurnFreq_leaveCluster,'Normalization','pdf','DisplayStyle','stairs')
        hold on
        histogram(omegaTurnFreq_loneWorm,'Normalization','pdf','DisplayStyle','stairs')
        legend(leaveClusterLegend, loneWormLegend)
        title([strains{strainCtr} '\_' wormnums{numCtr} '\_omegaTurnFreq'],'FontWeight','normal')
        xlabel('omega turn frequency (1/s)')
        ylabel('probability')
        xlim([0 0.8])
        ylim([0 18])
        set(omegaTurnFreqFig,'PaperUnits','centimeters')
        figurename = ['../figures/turns/omegaTurnFreq_' strains{strainCtr} '_' wormnums{numCtr} '_' phase '_CL'];
        savefig(omegaTurnFreqFig,[figurename '.fig'])
        %exportfig(omegaTurnFreqFig,[figurename '.eps'],exportOptions)
        %system(['epstopdf ' figurename '.eps']);
        %system(['rm ' figurename '.eps']);
    end
end