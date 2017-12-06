% since omega turn feature from the features.hdf5 is not very useful for
% detecting omega turns for example during leave cluster segments of the 
% trajectory (see comment on top of the analyzeTurnFeatures.m, this script 
% takes the individually isolated trajectories and
% runs the omega/upsilon turn detector on them separately. 

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
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_r_list.xlsx'],1,'A1:E15','basic');
        numFiles = length(filenames);
        % create cell arrays to hold individual movie values to be pooled
        %
        omegaTurnFreq_leaveCluster = cell(numFiles,1);
        omegaTurnFreq_loneWorm = cell(numFiles,1);
        omegaTurnFreqFig = figure;
        %
        omegaTurnFrameCount_all = cell(numFiles,1);
        omegaTurnFrameCount_leaveCluster = cell(numFiles,1);
        omegaTurnFrameCount_loneWorm = cell(numFiles,1);
        upsilonTurnFrameCount_all = cell(numFiles,1);
        upsilonTurnFrameCount_leaveCluster = cell(numFiles,1);
        upsilonTurnFrameCount_loneWorm = cell(numFiles,1);
        omegaTurnFrameCountFig = figure;
        upsilonTurnFrameCountFig = figure;
        %
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
            features.filtered = ismember(features.skeleton_id+1,find(trajData.filtered)); % use trajData.filtered to filter out unwanted data
            % find worms that have just left a cluster vs lone worms
            [leaveClusterLogInd, loneWormLogInd, ~,~] = findWormCategory(filename,inClusterNeighbourNum,minNeighbrDist,postExitDuration);
            leaveClusterLogIndFeat = ismember(features.skeleton_id+1,find(trajData.filtered & leaveClusterLogInd)); % make clusterClusterLogInd the same size as features.filtered
            loneWormLogIndFeat = ismember(features.skeleton_id+1,find(trajData.filtered & loneWormLogInd));
            
            %% calculate or extract desired feature values
            % detect omega turns by looping through each worm path
            %%%%%%%%%%%%%%%%%% really here we should break worm traj down
            %%%%%%%%%%%%%%%%%% into leave cluster and lone segments - see
            %%%%%%%%%%%%%%%%%% analyzeHeadTurn.m
            uniqueWormpaths = unique(trajData.worm_index_joined);
            worm_xcoords = squeeze(skelData(1,:,:))'; 
            worm_ycoords = squeeze(skelData(2,:,:))';
            omegaPaths_all = cell(numel(uniqueWormpaths),1);
            omegaPaths_leaveCluster = cell(numel(uniqueWormpaths),1);
            omegaPaths_loneWorm = cell(numel(uniqueWormpaths),1);
            upsilonPaths_all = cell(numel(uniqueWormpaths),1);
            upsilonPaths_leaveCluster = cell(numel(uniqueWormpaths),1);
            upsilonPaths_loneWorm = cell(numel(uniqueWormpaths),1);
            for wormpathCtr = 1:numel(uniqueWormpaths)
                wormpathLogInd = trajData.worm_index_joined==uniqueWormpaths(wormpathCtr);
                wormpath_xcoords_all = worm_xcoords(wormpathLogInd,:);
                wormpath_ycoords_all = worm_ycoords(wormpathLogInd,:);
                wormpath_xcoords_leaveCluster = worm_xcoords((wormpathLogInd & leaveClusterLogInd & trajData.filtered),:);
                wormpath_ycoords_leaveCluster = worm_ycoords((wormpathLogInd & leaveClusterLogInd & trajData.filtered),:);
                wormpath_xcoords_loneWorm = worm_xcoords(wormpathLogInd & loneWormLogInd & trajData.filtered,:);
                wormpath_ycoords_loneWorm = worm_ycoords(wormpathLogInd & loneWormLogInd & trajData.filtered,:);
                % randomly sample a section of the specified postExitDuration for trajectories
                if size(wormpath_xcoords_all,1)>(postExitDuration+1)*frameRate 
                    % postExitDuration + 1 to account for the extra second needed for smoothing later
                    firstFrame = randi((size(wormpath_xcoords_all,1)-(postExitDuration+1)*frameRate),1);
                    lastFrame = firstFrame + (postExitDuration+1)*frameRate;
                    wormpath_xcoords_all = wormpath_xcoords_all(firstFrame:lastFrame,:);
                    wormpath_ycoords_all = wormpath_ycoords_all(firstFrame:lastFrame,:);
                end
                if size(wormpath_xcoords_leaveCluster,1)>(postExitDuration+1)*frameRate
                    firstFrame = 1;
                    lastFrame = firstFrame + (postExitDuration+1)*frameRate;
                    wormpath_xcoords_leaveCluster = wormpath_xcoords_leaveCluster(firstFrame:lastFrame,:);
                    wormpath_ycoords_leaveCluster = wormpath_ycoords_leaveCluster(firstFrame:lastFrame,:);
                end
                if size(wormpath_xcoords_loneWorm,1)>(postExitDuration+1)*frameRate
                    firstFrame = randi(size(wormpath_xcoords_loneWorm,1)-(postExitDuration+1)*frameRate,1);
                    lastFrame = firstFrame + (postExitDuration+1)*frameRate;
                    wormpath_xcoords_loneWorm = wormpath_xcoords_loneWorm(firstFrame:lastFrame,:);
                    wormpath_ycoords_loneWorm = wormpath_ycoords_loneWorm(firstFrame:lastFrame,:);
                end
                [angleArray_all,meanAngles_all] = makeAngleArray(wormpath_xcoords_all,wormpath_ycoords_all);
                [angleArray_leaveCluster,meanAngles_leaveCluster] = makeAngleArray(wormpath_xcoords_leaveCluster,wormpath_ycoords_leaveCluster);
                [angleArray_loneWorm,meanAngles_loneWorm] = makeAngleArray(wormpath_xcoords_loneWorm,wormpath_ycoords_loneWorm);
                angleArray_all = angleArray_all + meanAngles_all;
                angleArray_leaveCluster = angleArray_leaveCluster + meanAngles_leaveCluster;
                angleArray_loneWorm = angleArray_loneWorm + meanAngles_loneWorm;
                headAngle_all = nanmean(angleArray_all(:,1:8),2); %restrict to head angles only
                headAngle_leaveCluster = nanmean(angleArray_leaveCluster(:,1:8),2);
                headAngle_loneWorm = nanmean(angleArray_loneWorm(:,1:8),2);
                % detect omega turns
                [omegaFrames_all,upsilonFrames_all] = omegaUpsilonDetectCurvature...
                    (angleArray_all',logical(zeros(size(angleArray_all,2),size(angleArray_all,1))));
                [omegaFrames_leaveCluster,upsilonFrames_leaveCluster] = omegaUpsilonDetectCurvature...
                    (angleArray_leaveCluster',logical(zeros(size(angleArray_leaveCluster,2),size(angleArray_leaveCluster,1))));
                [omegaFrames_loneWorm,upsilonFrames_loneWorm] = omegaUpsilonDetectCurvature...
                    (angleArray_loneWorm',logical(zeros(size(angleArray_loneWorm,2),size(angleArray_loneWorm,1))));
                omegaPaths_all{wormpathCtr}=omegaFrames_all(omegaFrames_all~=0);
                omegaPaths_leaveCluster{wormpathCtr}=omegaFrames_leaveCluster(omegaFrames_leaveCluster~=0);
                omegaPaths_loneWorm{wormpathCtr}=omegaFrames_loneWorm(omegaFrames_loneWorm~=0);
                upsilonPaths_all{wormpathCtr}=upsilonFrames_all(upsilonFrames_all~=0);
                upsilonPaths_leaveCluster{wormpathCtr}=upsilonFrames_leaveCluster(upsilonFrames_leaveCluster~=0);
                upsilonPaths_loneWorm{wormpathCtr}=upsilonFrames_loneWorm(upsilonFrames_loneWorm~=0);
            end
            %%%%%%%%%%%%%%%%% use full trajectory frequency %%%%%%%%%%%%%%%%%%%%%%%
            % write additional turn features into cell array for pooling across movies later
            uniqueLeaveClusterWorm = unique(features.worm_index(features.filtered&leaveClusterLogIndFeat));
            uniqueLoneWorm = unique(features.worm_index(features.filtered&loneWormLogIndFeat));
            omegaTurnFreq_leaveCluster_thisFile = zeros(1,numel(uniqueLeaveClusterWorm));
            omegaTurnFreq_loneWorm_thisFile = zeros(1,numel(uniqueLoneWorm));
            for wormpathCtr = 1:numel(uniqueLeaveClusterWorm)
                omegaTurnFreq_leaveCluster_thisFile(wormpathCtr)= h5read(strrep(filename,'skeletons','features'),['/features_events/worm_' num2str(uniqueLeaveClusterWorm(wormpathCtr)) '/omega_turns_frequency']);
            end
            for wormpathCtr = 1:numel(uniqueLoneWorm)
                omegaTurnFreq_loneWorm_thisFile(wormpathCtr)= h5read(strrep(filename,'skeletons','features'),['/features_events/worm_' num2str(uniqueLoneWorm(wormpathCtr)) '/omega_turns_frequency']);
            end
            omegaTurnFreq_leaveCluster{fileCtr} = omegaTurnFreq_leaveCluster_thisFile;
            omegaTurnFreq_loneWorm{fileCtr} = omegaTurnFreq_loneWorm_thisFile;
            %%%%%%%%%%%%%%%% use omega turn detector %%%%%%%%%%%%%%%%%%%%%
            % pool data from all wormpaths from a single file
            omegaPaths_all = vertcat(omegaPaths_all{:});
            omegaPaths_leaveCluster = vertcat(omegaPaths_leaveCluster{:});
            omegaPaths_loneWorm = vertcat(omegaPaths_loneWorm{:});
            upsilonPaths_all = vertcat(upsilonPaths_all{:});
            upsilonPaths_leaveCluster = vertcat(upsilonPaths_leaveCluster{:});
            upsilonPaths_loneWorm = vertcat(upsilonPaths_loneWorm{:});
            omegaTurnFrameCount_all{fileCtr} = omegaPaths_all;
            omegaTurnFrameCount_leaveCluster{fileCtr} = omegaPaths_leaveCluster;
            omegaTurnFrameCount_loneWorm{fileCtr} = omegaPaths_loneWorm;
            upsilonTurnFrameCount_all{fileCtr} = upsilonPaths_all;
            upsilonTurnFrameCount_leaveCluster{fileCtr} = upsilonPaths_leaveCluster;
            upsilonTurnFrameCount_loneWorm{fileCtr} = upsilonPaths_loneWorm;
            % update counter for final n number across movies
            leaveClusterWormCount = leaveClusterWormCount + numel(unique(features.worm_index(leaveClusterLogIndFeat)));
            loneWormCount = loneWormCount + numel(unique(features.worm_index(loneWormLogIndFeat)));
        end
        % pool data from all files belonging to the same strain and worm density
        omegaTurnFreq_leaveCluster = horzcat(omegaTurnFreq_leaveCluster{:})';
        omegaTurnFreq_loneWorm = horzcat(omegaTurnFreq_loneWorm{:})';
        %
        omegaTurnFrameCount_all = horzcat(omegaTurnFrameCount_all{:});
        omegaTurnFrameCount_leaveCluster = horzcat(omegaTurnFrameCount_leaveCluster{:});
        omegaTurnFrameCount_loneWorm = horzcat(omegaTurnFrameCount_loneWorm{:});
        upsilonTurnFrameCount_all = horzcat(upsilonTurnFrameCount_all{:});
        upsilonTurnFrameCount_leaveCluster = horzcat(upsilonTurnFrameCount_leaveCluster{:});
        upsilonTurnFrameCount_loneWorm = horzcat(upsilonTurnFrameCount_loneWorm{:});
        
        %% plot data, format, and export
        %omegaTurnFreq figure
        set(0,'CurrentFigure',omegaTurnFreqFig)
        histogram(omegaTurnFreq_leaveCluster,'Normalization','pdf','DisplayStyle','stairs')
        hold on
        histogram(omegaTurnFreq_loneWorm,'Normalization','pdf','DisplayStyle','stairs')
        legend('leave cluster', 'lone worm')
        title([strains{strainCtr} '\_' wormnums{numCtr} '\_omegaTurnFreq'],'FontWeight','normal')
        xlabel('omega turn frequency (1/s)')
        ylabel('probability')
        xlim([0 0.8])
        ylim([0 18])
        set(omegaTurnFreqFig,'PaperUnits','centimeters')
        figurename = ['figures/turns/omegaTurnFreq_' strains{strainCtr} '_' wormnums{numCtr} '_' phase '_CL'];
        %savefig(omegaTurnFreqFig,[figurename '.fig'])
        %exportfig(omegaTurnFreqFig,[figurename '.eps'],exportOptions)
        %system(['epstopdf ' figurename '.eps']);
        %system(['rm ' figurename '.eps']);
        
        %omegaTurnFrameDistFig
        set(0,'CurrentFigure',omegaTurnFrameCountFig)
        histogram(omegaTurnFrameCount_all,'Normalization','pdf','DisplayStyle','stairs')
        hold on
        histogram(omegaTurnFrameCount_leaveCluster,'Normalization','pdf','DisplayStyle','stairs')
        histogram(omegaTurnFrameCount_loneWorm,'Normalization','pdf','DisplayStyle','stairs')
        legend('all', 'leave cluster', 'lone worm')
        title([strains{strainCtr} '\_' wormnums{numCtr} '\_omegaTurnFrameDist'],'FontWeight','normal')
        set(omegaTurnFrameCountFig,'PaperUnits','centimeters')
        figurename = ['figures/turns/omegaTurnFrameCount_' strains{strainCtr} '_' wormnums{numCtr} '_' phase '_ACL'];
        
        %upsilonTurnFrameDistFig
        set(0,'CurrentFigure',upsilonTurnFrameCountFig)
        histogram(upsilonTurnFrameCount_all,'Normalization','pdf','DisplayStyle','stairs')
        hold on
        histogram(upsilonTurnFrameCount_leaveCluster,'Normalization','pdf','DisplayStyle','stairs')
        histogram(upsilonTurnFrameCount_loneWorm,'Normalization','pdf','DisplayStyle','stairs')
        legend('all', 'leave cluster', 'lone worm')
        title([strains{strainCtr} '\_' wormnums{numCtr} '\_upsilonTurnFrameDist'],'FontWeight','normal')
        set(omegaTurnFrameCountFig,'PaperUnits','centimeters')
        figurename = ['figures/turns/upsilonTurnFrameCount_' strains{strainCtr} '_' wormnums{numCtr} '_' phase '_ACL'];
    end
end