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

pixelsize = 100/19.5; % 100 microns are 19.5 pixels

strains = {'N2','npr1'};
wormnums = {'40'}%{'1W','40','HD'};
intensityThresholds_g = containers.Map({'40','HD','1W'},{60, 40, 100});
maxBlobSize_r = 2.5e5;
maxBlobSize_g = 1e4;
minNeighbrDist = 2000;
plotColors = lines(length(wormnums));

for strainCtr = 1:length(strains)
    revFreqFig = figure; hold on
    for numCtr = 1:length(wormnums)
        revDurFig = figure; hold on
        revInterTimeFig = figure; hold on
        wormnum = wormnums{numCtr};
        %% load data
        filenames_g = importdata(['datalists/' strains{strainCtr} '_' wormnum '_g_list.txt']);
        numFiles = length(filenames_g);
        reversalfreq_lone = NaN(numFiles,1);
        reversaldurations_lone = cell(numFiles,1);
        reversalfreq_incluster = NaN(numFiles,1);
        reversaldurations_incluster = cell(numFiles,1);
        reversalfreq_smallcluster = NaN(numFiles,1);
        reversaldurations_smallcluster = cell(numFiles,1);
        interrevT_lone = cell(numFiles,1);
        interrevT_incluster = cell(numFiles,1);
        interrevT_smallCluster = cell(numFiles,1);
        interrevT_lone_censored = cell(numFiles,1);
        interrevT_incluster_censored = cell(numFiles,1);
        interrevT_smallCluster_censored = cell(numFiles,1);
        frameRateAll = double(h5readatt(filenames_g{1},'/plate_worms','expected_fps')); % load one frameRate for use outside parfor loop
        parfor fileCtr = 1:numFiles % can be parfor
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
                min_neighbr_dist = h5read(filename_g,'/min_neighbr_dist');
                num_close_neighbrs = h5read(filename_g,'/num_close_neighbrs');
                neighbr_dist = h5read(filename_g,'/neighbr_distances');
                loneWormLogInd = min_neighbr_dist>=minNeighbrDist;
                inClusterLogInd = num_close_neighbrs>=3;
                smallClusterLogInd = (num_close_neighbrs==1 & neighbr_dist(:,2)>=minNeighbrDist)...
                    |(num_close_neighbrs==2 & neighbr_dist(:,3)>=minNeighbrDist)...
                    |(num_close_neighbrs==3 & neighbr_dist(:,4)>=minNeighbrDist)...
                    |(num_close_neighbrs==4 & neighbr_dist(:,5)>=minNeighbrDist);
            else
                loneWormLogInd = true(size(trajData_g.frame_number));
                inClusterLogInd = false(size(trajData_g.frame_number));
                smallClusterLogInd = false(size(trajData_g.frame_number));
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
            
            [ inclusterReversalsLogInd, interRevTimesCluster, revDurationCluster, interrevT_incluster_censored{fileCtr} ] = ...
                filterReversalsByClusterStatus(revStartInd, inClusterLogInd,...
                interRevTime, revDuration, incompleteInterRev);
            
            [ smallClusterReversalsLogInd, interRevTimesSmallCluster, revDurationSmallCluster, interrevT_smallCluster_censored{fileCtr} ] = ...
                filterReversalsByClusterStatus(revStartInd, smallClusterLogInd,...
                interRevTime, revDuration, incompleteInterRev);
            % subtracting revDuration will more accurately reflect the
            % interreversal time
            interrevT_lone{fileCtr} = (interRevTimesLone - revDurationLone)/frameRate;
            % ignore negative interRevTimes, as this (most likely) means
            % that the track was lost during a reversal
            interrevT_lone{fileCtr}(interrevT_lone{fileCtr}<0) = NaN;
            if ~strcmp(wormnum,'1W')
                interrevT_incluster{fileCtr} = (interRevTimesCluster - revDurationCluster)/frameRate;
                interrevT_incluster{fileCtr}(interrevT_incluster{fileCtr}<0) = NaN;
                interrevT_smallCluster{fileCtr} = (interRevTimesSmallCluster - revDurationSmallCluster)/frameRate;
                interrevT_smallCluster{fileCtr}(interrevT_smallCluster{fileCtr}<0) = NaN;
            end
            % counting reversal events
            reversalfreq_lone(fileCtr) = countReversalFrequency(loneReversalsLogInd,...
                frameRate, speedSigned, loneWormLogInd );
            reversalfreq_incluster(fileCtr) = countReversalFrequency(inclusterReversalsLogInd,...
                frameRate, speedSigned, inClusterLogInd );
            reversalfreq_smallcluster(fileCtr) =countReversalFrequency(smallClusterReversalsLogInd,...
                frameRate, speedSigned, smallClusterLogInd );
            reversaldurations_lone{fileCtr} = revDurationLone(~untrackedRevEnds(loneReversalsLogInd))/frameRate;
            reversaldurations_incluster{fileCtr} = revDurationCluster(~untrackedRevEnds(inclusterReversalsLogInd))/frameRate;
            reversaldurations_smallcluster{fileCtr} = revDurationSmallCluster(~untrackedRevEnds(smallClusterReversalsLogInd))/frameRate;
        end
        %pool data from all files
        interrevT_lone = vertcat(interrevT_lone{:});
        interrevT_incluster = vertcat(interrevT_incluster{:});
        interrevT_smallCluster = vertcat(interrevT_smallCluster{:});
        interrevT_lone_censored = vertcat(interrevT_lone_censored{:});
        interrevT_incluster_censored = vertcat(interrevT_incluster_censored{:});
        interrevT_smallCluster_censored = vertcat(interrevT_smallCluster_censored{:});
        %% plot data
        % inter-reversal time
        set(0,'CurrentFigure',revInterTimeFig)
        ecdf(interrevT_lone,'Bounds','on','function','survivor','censoring',interrevT_lone_censored)
        hold on
        if ~strcmp(wormnum,'1W')
%             ecdf(interrevT_smallCluster,'Bounds','on','function','survivor','censoring',interrevT_smallCluster_censored)
            ecdf(interrevT_incluster,'Bounds','on','function','survivor','censoring',interrevT_incluster_censored)
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
%             legend(revInterTimeFig.Children.Children([9 6 3]),{'lone worms','small cluster','in cluster'})
                        legend(revInterTimeFig.Children.Children([6 3]),{'lone worms','in cluster'})
        else
            legend(revInterTimeFig.Children,'single worms')
        end
        figurename = ['figures/reversals/reversalintertime_pharynx_' strains{strainCtr} '_' wormnum];
        exportfig(revInterTimeFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        % reversal frequency from counts
        set(0,'CurrentFigure',revFreqFig)
        notBoxPlot([reversalfreq_lone,reversalfreq_smallcluster,reversalfreq_incluster],...
            numCtr+[-0.3 0 0.3],'markMedian',true,'jitter',0.2)%,'style','line')
        
        revFreqFig.Children.XLim = [0 length(wormnums)+1];
        revFreqFig.Children.YLim = [0 0.55];
        % reversal durations
        reversaldurations_lone = vertcat(reversaldurations_lone{:});
        reversaldurations_incluster = vertcat(reversaldurations_incluster{:});
        reversaldurations_smallcluster = vertcat(reversaldurations_smallcluster{:});
        if strcmp(wormnum,'1W')
            histogram(revDurFig.Children,reversaldurations_lone,0:1/frameRateAll:15,...
                'Normalization','probability','DisplayStyle','stairs');
        else
            histogram(revDurFig.Children,reversaldurations_lone,0:3/frameRateAll:15,...
                'Normalization','probability','DisplayStyle','stairs');
%             histogram(revDurFig.Children,reversaldurations_smallcluster,0:3/frameRateAll:15,...
%                 'Normalization','probability','DisplayStyle','stairs','EdgeColor',0.5*ones(1,3));
            histogram(revDurFig.Children,reversaldurations_incluster,0:3/frameRateAll:15,...
                'Normalization','probability','DisplayStyle','stairs','EdgeColor','r');
        end
        revDurFig.Children.YTick = 0:0.1:0.5;
        title(revDurFig.Children,[strains{strainCtr} ' ' wormnum],'FontWeight','normal');
        set(revDurFig,'PaperUnits','centimeters')
        xlabel(revDurFig.Children,'time (s)')
        ylabel(revDurFig.Children,'P')
        revDurFig.Children.XLim = [0 8];
        if ~strcmp(wormnum,'1W')
%             legend(revDurFig.Children,{'lone worms','small cluster','in cluster'})
                        legend(revDurFig.Children,{'lone worms','in cluster'})
        else
            legend(revDurFig.Children,'single worms')
        end
        figurename = ['figures/reversals/reversaldurations_pharynx_' strains{strainCtr} '_' wormnum];
        exportfig(revDurFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
    end
    %% format and export figures
    title(revFreqFig.Children,strains{strainCtr},'FontWeight','normal');
    set(revFreqFig,'PaperUnits','centimeters')
    revFreqFig.Children.XTick = 1:length(wormnums);
    revFreqFig.Children.XTickLabel = strrep(wormnums,'HD','200');
    revFreqFig.Children.XLabel.String = 'worm number';
    revFreqFig.Children.YLabel.String = 'reversals (1/s)';
    revFreqFig.Children.YLim(1) = 0;
    figurename = ['figures/reversals/reversalfrequency_pharynx_' strains{strainCtr}];
    exportfig(revFreqFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end