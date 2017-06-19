% calculate various single worm statistics for different numbers of worms
% on plate

% issues / todo:

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
wormnums = {'1W','40','HD'};
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
        for fileCtr = 1:numFiles % can be parfor?
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
            %% calculate stats
            if ~strcmp(wormnum,'1W')
                min_neighbr_dist = h5read(filename_g,'/min_neighbr_dist');
                num_close_neighbrs = h5read(filename_g,'/num_close_neighbrs');
                neighbr_dist = h5read(filename_g,'/neighbr_distances');
                loneWormLogInd = min_neighbr_dist>=minNeighbrDist;
                loneWormInd = find(loneWormLogInd);
                inClusterLogInd = num_close_neighbrs>=3;
                inClusterInd = find(inClusterLogInd);
                smallClusterLogInd = (num_close_neighbrs==1 & neighbr_dist(:,2)>=minNeighbrDist)...
                    |(num_close_neighbrs==2 & neighbr_dist(:,3)>=minNeighbrDist)...
                    |(num_close_neighbrs==3 & neighbr_dist(:,4)>=minNeighbrDist)...
                    |(num_close_neighbrs==4 & neighbr_dist(:,5)>=minNeighbrDist);
                smallClusterInd = find(smallClusterLogInd);
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
            speedSigned(trajData_g.has_skeleton~=1)=NaN;
            % ignore skeletons otherwise filtered out
            speedSigned(~trajData_g.filtered) = NaN;
            % smooth speed to denoise
            speedSigned = smooth(speedSigned,3,'moving');
            % find reversals in midbody speed
            [revStartInd, revDuration, untrackedEnds] = findReversals(...
                speedSigned,trajData_g.worm_index_joined);
            % detect cluster status of reversals and features
            loneReversalsLogInd = ismember(revStartInd,find(loneWormLogInd));
            inclusterReversalsLogInd = ismember(revStartInd,find(inClusterLogInd));
            smallClusterReversalsLogInd = ismember(revStartInd,find(smallClusterLogInd));
            % estimate reversal frequencies and durations
            interRevTimesLone = diff(revStartInd(loneReversalsLogInd));
            interRevTimesCluster = diff(revStartInd(inclusterReversalsLogInd));
            interRevTimesSmallCluster = diff(revStartInd(smallClusterReversalsLogInd));
            revDurationLone = revDuration(loneReversalsLogInd);
            revDurationCluster = revDuration(inclusterReversalsLogInd);
            revDurationSmallCluster = revDuration(smallClusterReversalsLogInd);
            % censor those inter-reversal times that arise from
            % non-contiguous reversal sequences, eg when cluster status
            % changes btw reversals
            loneReversalInd = find(loneReversalsLogInd);
            nonContRevsLone = find(diff(loneReversalInd)~=1);
            for nCRctr = nonContRevsLone'
                % assign last time of same cluster status after each
                % non-contiguous reversal, and mark rev-intertime as censored
                interRevTimesLone(nCRctr) = ...
                    loneWormInd(find(loneWormInd<revStartInd(loneReversalInd(nCRctr)+1),1,'last'))...
                    -revStartInd(loneReversalInd(nCRctr));
            end
            interrevT_lone_censored{fileCtr} = false(size(interRevTimesLone));
            interrevT_lone_censored{fileCtr}(nonContRevsLone) = true;
            %
            inclusterReversalInd = find(inclusterReversalsLogInd);
            nonContRevsCluster = find(diff(inclusterReversalInd)~=1);
            for nCRctr = nonContRevsCluster'
                interRevTimesCluster(nCRctr) = ...
                    inclusterReversalInd(find(inclusterReversalInd<revStartInd(inclusterReversalInd(nCRctr)+1),1,'last'))...
                    -revStartInd(inclusterReversalInd(nCRctr));
            end
            interrevT_incluster_censored{fileCtr} = false(size(interRevTimesCluster));
            interrevT_incluster_censored{fileCtr}(nonContRevsCluster) = true;
            %
            smallClusterReversalInd = find(smallClusterReversalsLogInd);
            nonContRevsSmallCluster = find(diff(smallClusterReversalInd)~=1);
            for nCRctr = nonContRevsSmallCluster'
                interRevTimesSmallCluster(nCRctr) = ...
                    smallClusterReversalInd(find(smallClusterReversalInd<revStartInd(smallClusterReversalInd(nCRctr)+1),1,'last'))...
                    -revStartInd(smallClusterReversalInd(nCRctr));
            end
            interrevT_smallCluster_censored{fileCtr} = false(size(interRevTimesSmallCluster));
            interrevT_smallCluster_censored{fileCtr}(nonContRevsSmallCluster) = true;
%             % censor those inter-reversal times where the worm_index
%             % changes
%             wormChangeLogInd = trajData_g.worm_index_joined(revStartInd(loneReversalInd(1:end-1)))~=...
%                 trajData_g.worm_index_joined(revStartInd(loneReversalInd(1:end-1))+interRevTimesLone);
            
            % subtracting revDuration will more accurately reflect the
            % inter reversal time
            interrevT_lone{fileCtr} = (interRevTimesLone - revDurationLone(1:end-1))/frameRate;
             % ignore negative interRevTimes, as this (most likely) means
            % that the track was lost during a reversal
            interrevT_lone{fileCtr}(interrevT_lone{fileCtr}<0) = NaN;
            interrevT_lone_censored{fileCtr} = interrevT_lone_censored{fileCtr}&...
                untrackedEnds(find(loneReversalsLogInd,numel(interRevTimesLone),'first'));
            if ~strcmp(wormnum,'1W')
                interrevT_incluster{fileCtr} = (interRevTimesCluster - revDurationCluster(1:end-1))/frameRate;
                            interrevT_incluster{fileCtr}(interrevT_incluster{fileCtr}<0) = NaN;
                interrevT_smallCluster{fileCtr} = (interRevTimesSmallCluster - revDurationSmallCluster(1:end-1))/frameRate;
                            interrevT_smallCluster{fileCtr}(interrevT_smallCluster{fileCtr}<0) = NaN;
                interrevT_incluster_censored{fileCtr} = interrevT_incluster_censored{fileCtr}&...
                    untrackedEnds(find(inclusterReversalsLogInd,numel(interRevTimesCluster),'first'));
                interrevT_smallCluster_censored{fileCtr} = interrevT_smallCluster_censored{fileCtr}&...
                    untrackedEnds(find(smallClusterReversalsLogInd,numel(interRevTimesSmallCluster),'first'));
            end
            % counting reversal events
            Nrev_lone = nnz(loneReversalsLogInd);
            Nrev_incluster = nnz(inclusterReversalsLogInd);
            Nrev_smallCluster = nnz(smallClusterReversalsLogInd);
            T_lone = nnz(loneWormLogInd)/frameRate;
            T_incluster = nnz(inClusterLogInd)/frameRate;
            T_smallcluster = nnz(smallClusterLogInd)/frameRate;
            Trev_lone = nnz(speedSigned(loneWormLogInd)<0)/frameRate;
            Trev_incluster = nnz(speedSigned(inClusterLogInd)<0)/frameRate;
            Trev_smallcluster = nnz(speedSigned(smallClusterReversalsLogInd)<0)/frameRate;
            reversalfreq_lone(fileCtr) = Nrev_lone./(T_lone - Trev_lone);
            reversalfreq_incluster(fileCtr) = Nrev_incluster./(T_incluster - Trev_incluster);
            reversalfreq_smallcluster(fileCtr) = Nrev_smallCluster./(T_smallcluster - Trev_smallcluster);
            reversaldurations_lone{fileCtr} = revDurationLone(~untrackedEnds(loneReversalsLogInd))/frameRate;
            reversaldurations_incluster{fileCtr} = revDurationCluster(~untrackedEnds(inclusterReversalsLogInd))/frameRate;
            reversaldurations_smallcluster{fileCtr} = revDurationSmallCluster(~untrackedEnds(smallClusterReversalsLogInd))/frameRate;
        end
        %pool data from all files
        interrevT_lone = vertcat(interrevT_lone{:});
        interrevT_incluster = vertcat(interrevT_incluster{:});
        interrevT_smallCluster = vertcat(interrevT_smallCluster{:});
        interrevT_lone_censored = vertcat(interrevT_lone_censored{:});
        interrevT_incluster_censored = vertcat(interrevT_incluster_censored{:});
        interrevT_smallCluster_censored = vertcat(interrevT_smallCluster_censored{:});
        %% plot data
        set(0,'CurrentFigure',revInterTimeFig)
        ecdf(interrevT_lone,'Bounds','on','function','survivor','censoring',interrevT_lone_censored)
        hold on
        if ~strcmp(wormnum,'1W')
            ecdf(interrevT_smallCluster,'Bounds','on','function','survivor','censoring',interrevT_smallCluster_censored)
            ecdf(interrevT_incluster,'Bounds','on','function','survivor','censoring',interrevT_incluster_censored)
        end
        set(revInterTimeFig.Children,'YScale','log')
        title(revInterTimeFig.Children,[strains{strainCtr} ' ' wormnum],'FontWeight','normal');
        set(revInterTimeFig,'PaperUnits','centimeters')
        revInterTimeFig.Children.XLabel.String = 'inter-reversal time (s)';
        revInterTimeFig.Children.YLabel.String = 'cumulative probability';
        revInterTimeFig.Children.XLim(2) = 60;
        if ~strcmp(wormnum,'1W')
            revInterTimeFig.Children.YLim(1) = 1e-4;
            legend(revInterTimeFig.Children.Children([9 6 3]),{'lone worms','small cluster','in cluster'})
            %             legend(revInterTimeFig.Children.Children([6 3]),{'lone worms','in cluster'})
        else
            legend(revInterTimeFig.Children,'single worms')
        end
        figurename = ['figures/reversals/reversalintertime_pharynx_' strains{strainCtr} '_' wormnum];
        exportfig(revInterTimeFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        %
        set(0,'CurrentFigure',revFreqFig)
        notBoxPlot([reversalfreq_lone,reversalfreq_smallcluster,reversalfreq_incluster],...
            numCtr+[-0.3 0 0.3],'markMedian',true,'jitter',0.2)%,'style','line')
        
        revFreqFig.Children.XLim = [0 length(wormnums)+1];
        revFreqFig.Children.YLim = [0 0.55];
        %
        reversaldurations_lone = vertcat(reversaldurations_lone{:});
        reversaldurations_incluster = vertcat(reversaldurations_incluster{:});
        reversaldurations_smallcluster = vertcat(reversaldurations_smallcluster{:});
        if strcmp(wormnum,'1W')
            histogram(revDurFig.Children,reversaldurations_lone,0:1/frameRate:15,...
                'Normalization','probability','DisplayStyle','stairs');
        else
            histogram(revDurFig.Children,reversaldurations_lone,0:3/frameRate:15,...
                'Normalization','probability','DisplayStyle','stairs');
            histogram(revDurFig.Children,reversaldurations_smallcluster,0:3/frameRate:15,...
                'Normalization','probability','DisplayStyle','stairs','EdgeColor',0.5*ones(1,3));
            histogram(revDurFig.Children,reversaldurations_incluster,0:3/frameRate:15,...
                'Normalization','probability','DisplayStyle','stairs','EdgeColor','r');
        end
        revDurFig.Children.YTick = 0:0.1:0.5;
        %
        title(revDurFig.Children,[strains{strainCtr} ' ' wormnum],'FontWeight','normal');
        set(revDurFig,'PaperUnits','centimeters')
        xlabel(revDurFig.Children,'time (s)')
        ylabel(revDurFig.Children,'P')
        revDurFig.Children.XLim = [0 10];
        if ~strcmp(wormnum,'1W')
            legend(revDurFig.Children,{'lone worms','small cluster','in cluster'})
            %             legend(revDurFig.Children,{'lone worms','in cluster'})
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