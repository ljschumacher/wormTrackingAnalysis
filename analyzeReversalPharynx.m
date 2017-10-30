function [] = analyzeReversalPharynx(dataset,phase,wormnum)
% calculate various single worm statistics for different numbers of worms
% on plate
% INPUTS
% dataset: 1 or 2. To specify which dataset to run the script for.
% phase: 'joining', 'fullMovie', or 'sweeping'. Script defines stationary phase as: starts at 10% into the movie, and stops at 60% into the movie (HA and N2) or at specified stopping frames (npr-1).
% wormnum: '1W', '40', or 'HD'
% OUTPUTS
% none returned, but figures are exported

% issues / todo:
% - for the kaplan-meier survival curves, we should not plot the inter rev
% time for lone worms, but the time from a random starting point until rev


exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

%% set parameters
if dataset ==1
    strains = {'npr1','N2'}; %{'npr1','HA','N2'}
elseif dataset ==2
    strains = {'npr1','N2'}; %{'npr1','N2'}
end
if dataset == 1
    intensityThresholds_g = containers.Map({'40','HD','1W'},{50, 40, 100});
elseif dataset ==2
    intensityThresholds_g = containers.Map({'40','HD','1W'},{60, 40, 100});
end
maxBlobSize_g = 1e4;
minNeighbrDist = 2000;% in microns
minPathLength = 50; % minimum path length of reversals to be included
inClusterNeighbourNum = 3;
postExitDuration = 10; % set the duration (in seconds) after a worm exits a cluster to be included in the leave cluster analysis
pixelsize = 100/19.5; % 100 microns are 19.5 pixels

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    revFreqFig = figure; hold on
    revInterTimeFig = figure; hold on
    %% load data
    if dataset ==1
        [phaseFrames,filenames_g,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_list.xlsx'],1,'A1:E15','basic');
    elseif dataset ==2
        [phaseFrames,filenames_g,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_g_list.xlsx'],1,'A1:E15','basic');
    end
    numFiles = length(filenames_g);
    reversalfreq_lone = NaN(numFiles,1);
    reversalfreq_leaveCluster = NaN(numFiles,1);
    reversalfreq_revCluster = NaN(numFiles,1);
    interrevT_lone = cell(numFiles,1);
    timeToRev_leaveCluster = cell(numFiles,1);
    timeToFwd_revCluster = cell(numFiles,1);
    interrevT_lone_censored = cell(numFiles,1);
    timeToRev_leaveCluster_censored = cell(numFiles,1);
    timeToFwd_revCluster_censored = cell(numFiles,1);
    for fileCtr = 1:numFiles % can be parfor
        filename_g = filenames_g{fileCtr};
        trajData_g = h5read(filename_g,'/trajectories_data');
        blobFeats_g = h5read(filename_g,'/blob_features');
        skelData_g = h5read(filename_g,'/skeleton');
        assert(size(skelData_g,1)==2&&size(skelData_g,2)==2,['unexpected skeleton size for ' filename_g]);
        frameRate = double(h5readatt(filename_g,'/plate_worms','expected_fps'));
        if frameRate == 0
            warning(['frame rate is zero for ' filename_g])
        end
        % filter by blob size and intensity
        trajData_g.filtered = filterIntensityAndSize(blobFeats_g,pixelsize,...
            intensityThresholds_g(wormnum),maxBlobSize_g);
        trajData_g.has_skeleton = squeeze(~any(any(isnan(skelData_g)))); % reset skeleton flag for pharynx data
        % check worm-indices are monotonically increasing
        assert(~any(diff(trajData_g.worm_index_joined)<0),['worm indices are not sorted as expected for ' filename_g])
        %% calculate stats
        if ~strcmp(wormnum,'1W')
            % find leave cluster and lone worms
            [leaveClusterLogInd, loneWormLogInd,~,~] = findWormCategory(filename_g,inClusterNeighbourNum,minNeighbrDist,postExitDuration);
            % apply phase restriction
            [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr);
            phaseFrameLogInd = trajData_g.frame_number < lastFrame & trajData_g.frame_number > firstFrame;
            loneWormLogInd(~phaseFrameLogInd) = false;
            leaveClusterLogInd(~phaseFrameLogInd) = false;
        else
            loneWormLogInd = true(size(trajData_g.frame_number));
            leaveClusterLogInd = false(size(trajData_g.frame_number));
        end
        %% load signed speed from blobFeats
        if any(phaseFrameLogInd)
            % sign speed based on relative orientation of velocity to midbody
            speedSigned = blobFeats_g.signed_speed*pixelsize*frameRate;
            % ignore first and last frames of each worm's track
            wormChangeIndcs = gradient(double(trajData_g.worm_index_joined))~=0;
            speedSigned(wormChangeIndcs)=NaN;
            % ignore frames with bad skeletonization
            speedSigned(~trajData_g.has_skeleton)=NaN;
            % ignore skeletons otherwise filtered out
            speedSigned(~trajData_g.filtered) = NaN;
            % ignore frames outside of specified phase
            if ~strcmp(phase,'fullMovie')
                speedSigned(~phaseFrameLogInd) =NaN;
            end
            % smooth speed to denoise
            speedSigned = smoothdata(speedSigned,'movmean',3,'omitnan');
            % find reversals in signed speed
            [revStartInd, revDuration, untrackedRevEnds, interRevTime, incompleteInterRev] = ...
                findReversals(speedSigned,trajData_g.worm_index_joined,minPathLength,frameRate);
%             [fwdStartInd, ~, ~, ~, ~] = ...
%                 findReversals(-speedSigned,trajData_g.worm_index_joined,minPathLength,frameRate);
            fwdStartInd = find(speedSigned(1:end-1)<0&speedSigned(2:end)>0);
            % if we subtract rev duration from interrevtime (below), we
            % need to set all reversals with untracked ends as incomplete interrevs
            incompleteInterRev = incompleteInterRev|untrackedRevEnds;
            % detect cluster status of reversals and associated features
            [ loneReversalsLogInd, interRevTimesLone, revDurationLone, interrevT_lone_censored{fileCtr} ] = ...
                filterReversalsByClusterStatus(revStartInd, loneWormLogInd,...
                interRevTime, revDuration, incompleteInterRev);
            % filter reversals after leaving cluster
            [ timeToRevLeaveCluster, timeToRev_leaveCluster_censored{fileCtr} ] = ...
                filterReversalsByEvent(revStartInd, leaveClusterLogInd, trajData_g.worm_index_joined,postExitDuration*frameRate);
            % filter reversals after reversing out of cluster
            [ timeToFwdRevCluster, timeToFwd_revCluster_censored{fileCtr} ] = ...
                filterReversalsByEvent(fwdStartInd, leaveClusterLogInd&speedSigned<0, trajData_g.worm_index_joined,postExitDuration*frameRate);
            % subtracting revDuration will more accurately reflect the
            % interreversal time
            interrevT_lone{fileCtr} = (interRevTimesLone - revDurationLone)/frameRate;
            % ignore negative interRevTimes, as this (most likely) means
            % that the track was lost during a reversal
            interrevT_lone{fileCtr}(interrevT_lone{fileCtr}<0) = NaN;
            if ~strcmp(wormnum,'1W')
                timeToRev_leaveCluster{fileCtr} = timeToRevLeaveCluster/frameRate;
                timeToFwd_revCluster{fileCtr} = timeToFwdRevCluster/frameRate;
            end
            % counting reversal events
            reversalfreq_lone(fileCtr) = countReversalFrequency(loneReversalsLogInd,...
                frameRate, speedSigned, loneWormLogInd );
            % count how many leave-cluster trajs end in reversals (are not
            % censored), and divide by the total time until reversals (or
            % trajectories end)
            reversalfreq_leaveCluster(fileCtr) = nnz(timeToRev_leaveCluster_censored{fileCtr})/sum(timeToRev_leaveCluster{fileCtr});
            reversalfreq_revCluster(fileCtr) = nnz(timeToFwd_revCluster_censored{fileCtr})/sum(timeToFwd_revCluster{fileCtr});
        end
    end
    %pool data from all files
    interrevT_lone = vertcat(interrevT_lone{:});
    timeToRev_leaveCluster = vertcat(timeToRev_leaveCluster{:});
    timeToFwd_revCluster = vertcat(timeToFwd_revCluster{:});
    interrevT_lone_censored = vertcat(interrevT_lone_censored{:});
    timeToRev_leaveCluster_censored = vertcat(timeToRev_leaveCluster_censored{:});
    timeToFwd_revCluster_censored = vertcat(timeToFwd_revCluster_censored{:});
    %% plot data, format and export figures
    % inter-reversal time
    set(0,'CurrentFigure',revInterTimeFig)
    ecdf(interrevT_lone,'Bounds','on','function','survivor','censoring',interrevT_lone_censored)
    hold on
    if ~strcmp(wormnum,'1W')
        ecdf(timeToRev_leaveCluster,'Bounds','on','function','survivor','censoring',timeToRev_leaveCluster_censored)
        ecdf(timeToFwd_revCluster,'Bounds','on','function','survivor','censoring',timeToFwd_revCluster_censored)
    end
    set(revInterTimeFig.Children,'YScale','log')
    title(revInterTimeFig.Children,[strains{strainCtr} ' ' wormnum],'FontWeight','normal');
    set(revInterTimeFig,'PaperUnits','centimeters')
    revInterTimeFig.Children.XLabel.String = 'inter-reversal time (s)';
    revInterTimeFig.Children.YLabel.String = 'cumulative probability';
    revInterTimeFig.Children.XLim(2) = 20;
    revInterTimeFig.Children.YLim(1) = 1e-1;
    if ~strcmp(wormnum,'1W')
        legend(revInterTimeFig.Children.Children([9 6 3]),{'lone worms','leaving cluster','rev-leaving cluster'})
    else
        legend(revInterTimeFig.Children,'single worms')
    end
    figurename = ['figures/reversals/phaseSpecific/reversalintertime_pharynx_'...
        strains{strainCtr} '_' wormnum '_' phase '_data' num2str(dataset)  '_jointraj' '_censored'];
    exportfig(revInterTimeFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
    
    % reversal frequency from counts
    set(0,'CurrentFigure',revFreqFig)
    notBoxPlot([reversalfreq_lone,reversalfreq_leaveCluster,reversalfreq_revCluster],...
        [-0.3 0 0.3],'markMedian',true,'jitter',0.2)%,'style','line')
    title(revFreqFig.Children,strains{strainCtr},'FontWeight','normal');
    set(revFreqFig,'PaperUnits','centimeters')
    if ~strcmp(wormnum,'1W')
        revFreqFig.Children.XTickLabel = {'lone','leave','rev-leave'};
    end
    revFreqFig.Children.XLabel.String = 'worm categories';
    revFreqFig.Children.YLabel.String = 'reversals (1/s)';
    revFreqFig.Children.YLim(1) = 0;
    revFreqFig.Children.YLim(2) = 1;
    figurename = ['figures/reversals/phaseSpecific/reversalfrequency_pharynx_'...
        strains{strainCtr} '_' wormnum '_' phase '_data' num2str(dataset) '_jointraj'];
    exportfig(revFreqFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end