% calculate various single worm statistics for different numbers of worms
% on plate

% issues / todo:
% - censoring reversals when a cluster status changes drastically reduces
% the length of observed interrev times

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
dataset = 2; % '1' or '2'. To specify which dataset to run the script for.
phase = 'stationary'; % 'fullMovie' or 'stationary'. Script defines stationary phase as: starts at 10% into the movie, and stops at 60% into the movie (HA and N2) or at specified stopping frames (npr-1).
if dataset ==1
    strains = {'npr1','N2'}; %{'npr1','HA','N2'}
elseif dataset ==2
    strains = {'npr1','N2'}; %{'npr1','N2'}
end
wormnums = {'40'}; %{'1W','40','HD'};
if dataset == 1
    intensityThresholds_g = containers.Map({'40','HD','1W'},{50, 40, 100});
elseif dataset ==2
    intensityThresholds_g = containers.Map({'40','HD','1W'},{60, 40, 100});
end
maxBlobSize_g = 1e4;
minNeighbrDist = 2000;
inClusterNeighbourNum = 3;
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
plotColors = lines(length(wormnums));

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    revFreqFig = figure; hold on
    for numCtr = 1:length(wormnums)
         revInterTimeFig = figure; hold on
        wormnum = wormnums{numCtr};
        %% load data
        if dataset ==1
            [lastFrames,filenames_g,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_list.xlsx'],1,'A1:B15','basic');
        elseif dataset ==2
            [lastFrames,filenames_g,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_g_list.xlsx'],1,'A1:B15','basic');
        end
        numFiles = length(filenames_g);
        reversalfreq_lone = NaN(numFiles,1);
        reversalfreq_notCluster = NaN(numFiles,1);
        interrevT_lone = cell(numFiles,1);
        interrevT_notCluster = cell(numFiles,1);
        interrevT_lone_censored = cell(numFiles,1);
        interrevT_notCluster_censored = cell(numFiles,1);
        frameRateAll = double(h5readatt(filenames_g{1},'/plate_worms','expected_fps')); % load one frameRate for use outside parfor loop
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
            if strcmp(phase,'fullMovie')
                lastFrame = double(max(trajData_g.frame_number));
            elseif strcmp(phase,'stationary')
                lastFrame = lastFrames(fileCtr);
            end
            % filter by blob size and intensity
            trajData_g.filtered = filterIntensityAndSize(blobFeats_g,pixelsize,...
                intensityThresholds_g(wormnum),maxBlobSize_g);
            trajData_g.has_skeleton = squeeze(~any(any(isnan(skelData_g)))); % reset skeleton flag for pharynx data
            % check worm-indices are monotonically increasing
            assert(~any(diff(trajData_g.worm_index_joined)<0),['worm indices are not sorted as expected for ' filename_g])
            %% calculate stats
            if ~strcmp(wormnum,'1W')
                min_neighbr_dist = h5read(filename_g,'/min_neighbr_dist');
                num_close_neighbrs = h5read(filename_g,'/num_close_neighbrs');
                neighbr_dist = h5read(filename_g,'/neighbr_distances');
                loneWormLogInd = min_neighbr_dist>=minNeighbrDist;
                inClusterLogInd = num_close_neighbrs>=inClusterNeighbourNum;
                % find frames when worms are leaving cluster
                leaveClusterLogInd = [false; inClusterLogInd(1:end-1)&~inClusterLogInd(2:end)];
                % add 5 more second after each leaving event
                leaveClusterExtendedInd = unique(find(leaveClusterLogInd) + [0:5*frameRate]);
                leaveClusterExtendedInd = leaveClusterExtendedInd(leaveClusterExtendedInd<numel(leaveClusterLogInd)); % exclude frames beyond highest frame number
                leaveClusterLogInd(leaveClusterExtendedInd) = true; 
                leaveClusterLogInd(inClusterLogInd) = false; % exclude times when worm moved back
                leaveClusterLogInd(loneWormLogInd)=false; % exclude worms that have become lone worm
                if strcmp(phase,'stationary')
                    % restrict movie to stationary phase
                    firstFrame = double(round(max(trajData_g.frame_number)/10)); % define starting frame as 10% into the movie
                    phaseFrameLogInd = trajData_g.frame_number < lastFrame & trajData_g.frame_number > firstFrame;
                    loneWormLogInd(~phaseFrameLogInd) = false;
                    leaveClusterLogInd(~phaseFrameLogInd) = false;
                end
            else
                loneWormLogInd = true(size(trajData_g.frame_number));
                leaveClusterLogInd = false(size(trajData_g.frame_number));
            end
            %% load signed speed from blobFeats
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
            if strcmp(phase,'stationary')
                speedSigned(~phaseFrameLogInd) =NaN;
            end
            % smooth speed to denoise
            speedSigned = smooth(speedSigned,3,'moving');
            % find reversals in midbody speed
            [revStartInd, revDuration, untrackedRevEnds, interRevTime, incompleteInterRev] = ...
                findReversals(speedSigned,trajData_g.worm_index_joined);
            % if we subtract rev duration from interrevtime (below), we
            % need to set all reversals with untracked ends as incomplete interrevs
            incompleteInterRev = incompleteInterRev|untrackedRevEnds;
            % detect cluster status of reversals and features
            [ loneReversalsLogInd, interRevTimesLone, revDurationLone, interrevT_lone_censored{fileCtr} ] = ...
                filterReversalsByClusterStatus(revStartInd, loneWormLogInd,...
                interRevTime, revDuration, incompleteInterRev);
            
            [ leaveClusterReversalsLogInd, interRevTimesCluster, revDurationCluster, interrevT_notCluster_censored{fileCtr} ] = ...
                filterReversalsByClusterStatus(revStartInd, leaveClusterLogInd,...
                interRevTime, revDuration, incompleteInterRev);
            % subtracting revDuration will more accurately reflect the
            % interreversal time
            interrevT_lone{fileCtr} = (interRevTimesLone - revDurationLone)/frameRate;
            % ignore negative interRevTimes, as this (most likely) means
            % that the track was lost during a reversal
            interrevT_lone{fileCtr}(interrevT_lone{fileCtr}<0) = NaN;
            if ~strcmp(wormnum,'1W')
                interrevT_notCluster{fileCtr} = (interRevTimesCluster - revDurationCluster)/frameRate;
                interrevT_notCluster{fileCtr}(interrevT_notCluster{fileCtr}<0) = NaN;
            end
            % counting reversal events
            reversalfreq_lone(fileCtr) = countReversalFrequency(loneReversalsLogInd,...
                frameRate, speedSigned, loneWormLogInd );
            reversalfreq_notCluster(fileCtr) = countReversalFrequency(leaveClusterReversalsLogInd,...
                frameRate, speedSigned, leaveClusterLogInd );
        end
        %pool data from all files
        interrevT_lone = vertcat(interrevT_lone{:});
        interrevT_notCluster = vertcat(interrevT_notCluster{:});
        interrevT_lone_censored = vertcat(interrevT_lone_censored{:});
        interrevT_notCluster_censored = vertcat(interrevT_notCluster_censored{:});
        %% plot data
        % inter-reversal time
        set(0,'CurrentFigure',revInterTimeFig)
        ecdf(interrevT_lone,'Bounds','on','function','survivor','censoring',interrevT_lone_censored)
        hold on
        if ~strcmp(wormnum,'1W')
            ecdf(interrevT_notCluster,'Bounds','on','function','survivor','censoring',interrevT_notCluster_censored)
        end
        set(revInterTimeFig.Children,'YScale','log')
        title(revInterTimeFig.Children,[strains{strainCtr} ' ' wormnum],'FontWeight','normal');
        set(revInterTimeFig,'PaperUnits','centimeters')
        revInterTimeFig.Children.XLabel.String = 'inter-reversal time (s)';
        revInterTimeFig.Children.YLabel.String = 'cumulative probability';
        %         revInterTimeFig.Children.XLim(2) = 30;
        %         revInterTimeFig.Children.YLim(1) = 0.1;
        revInterTimeFig.Children.XLim(2) = 60;
        revInterTimeFig.Children.YLim(1) = 1e-3;
        if ~strcmp(wormnum,'1W')
            legend(revInterTimeFig.Children.Children([6 3]),{'lone worms','not in cluster'})
        else
            legend(revInterTimeFig.Children,'single worms')
        end
        
        figurename = ['figures/reversals/phaseSpecific/reversalintertime_pharynx_' strains{strainCtr} '_' wormnum '_' phase '_data' num2str(dataset) '_censored'];
        exportfig(revInterTimeFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        % reversal frequency from counts
        set(0,'CurrentFigure',revFreqFig)
        notBoxPlot([reversalfreq_lone,reversalfreq_notCluster],...
            numCtr+[-0.2 0.2],'markMedian',true,'jitter',0.2)%,'style','line')
        
        revFreqFig.Children.XLim = [0 length(wormnums)+1];
        revFreqFig.Children.YLim = [0 0.55];
    end
    %% format and export figures
    title(revFreqFig.Children,strains{strainCtr},'FontWeight','normal');
    set(revFreqFig,'PaperUnits','centimeters')
    revFreqFig.Children.XTick = 1:length(wormnums);
    revFreqFig.Children.XTickLabel = strrep(wormnums,'HD','200');
    revFreqFig.Children.XLabel.String = 'worm number';
    revFreqFig.Children.YLabel.String = 'reversals (1/s)';
    revFreqFig.Children.YLim(1) = 0;
    figurename = ['figures/reversals/phaseSpecific/reversalfrequency_pharynx_' strains{strainCtr} '_' phase '_data' num2str(dataset) '_censored'];
    exportfig(revFreqFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end