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
minSkelLength_r = 850;
maxSkelLength_r = 1500;
minNeighbrDist = 1500;
midbodyIndcs = 19:33;
plotColors = lines(length(wormnums));
loneClusterRadius = 2000; % for defining small clusters

for strainCtr = 1:length(strains)
    revFreqFig = figure; hold on
    for numCtr = 1:length(wormnums)
        revDurFig = figure; hold on
        revInterTimeFig = figure; hold on
        wormnum = wormnums{numCtr};
        %% load data
        filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_r_list.txt']);
        if ~strcmp(wormnum,'1W')
            filenames_g = importdata(['datalists/' strains{strainCtr} '_' wormnum '_g_list.txt']);
        else
            filenames_g = {};
        end
        numFiles = length(filenames);
        reversalfreq_lone = NaN(numFiles,1);
        reversaldurations_lone = cell(numFiles,1);
        reversalfreq_inCluster = NaN(numFiles,1);
        reversaldurations_inCluster = cell(numFiles,1);
        reversalfreq_smallCluster = NaN(numFiles,1);
        reversaldurations_smallCluster = cell(numFiles,1);
        interrevT_lone = cell(numFiles,1);
        interrevT_inCluster = cell(numFiles,1);
        interrevT_smallCluster = cell(numFiles,1);
        for fileCtr = 1:numFiles % can be parfor?
            filename_r = filenames{fileCtr};
            trajData_r = h5read(filename_r,'/trajectories_data');
            blobFeats_r = h5read(filename_r,'/blob_features');
            skelData_r = h5read(filename_r,'/skeleton');
            % % featData = h5read(strrep(filename,'skeletons','features'),'/features_timeseries');
            frameRate = double(h5readatt(filename_r,'/plate_worms','expected_fps'));
            % filter red data by blob size and intensity
            if contains(filename_r,'55')||contains(filename_r,'54')
                intensityThreshold_r = 80;
            else
                intensityThreshold_r = 40;
            end
            trajData_r.filtered = filterIntensityAndSize(blobFeats_r,pixelsize,...
                intensityThreshold_r,maxBlobSize_r);
            % filter red data by skeleton length
            trajData_r.filtered = trajData_r.filtered&logical(trajData_r.is_good_skel)&...
                filterSkelLength(skelData_r,pixelsize,minSkelLength_r,maxSkelLength_r);
            %% calculate stats
            if ~strcmp(wormnum,'1W')
                min_neighbr_dist_rr = h5read(filename_r,'/min_neighbr_dist_rr');
                min_neighbr_dist_rg = h5read(filename_r,'/min_neighbr_dist_rg');
                num_close_neighbrs_rg = h5read(filename_r,'/num_close_neighbrs_rg')';
                neighbr_dist = h5read(filename_r,'/neighbr_distances');
                loneWorms = min_neighbr_dist_rr>=1100&min_neighbr_dist_rg>=1600;
                inCluster = num_close_neighbrs_rg>=3;
                smallCluster = trajData_r.filtered&...
                    ((num_close_neighbrs_rg == 2 & neighbr_dist(:,3)>=(loneClusterRadius))...
                    |(num_close_neighbrs_rg == 3 & neighbr_dist(:,4)>=(loneClusterRadius))...
                    |(num_close_neighbrs_rg == 4 & neighbr_dist(:,5)>=(loneClusterRadius)));
                % define lone, in cluster, and small cluster worms
            else
                loneWorms = true(size(trajData_r.frame_number));
                inCluster = false(size(trajData_r.frame_number));
                smallCluster = false(size(trajData_r.frame_number));
            end
            %% calculate overall midbody speeds, until we can link trajectories and features from the tracker
            % centroids of midbody skeleton
            midbody_x = mean(squeeze(skelData_r(1,midbodyIndcs,:)))*pixelsize;
            midbody_y = mean(squeeze(skelData_r(2,midbodyIndcs,:)))*pixelsize;
            % change in centroid position over time (issue of worm shifts?)
            dmidbody_xdt = gradient(midbody_x)*frameRate;
            dmidbody_ydt = gradient(midbody_y)*frameRate;
            % midbody speed and velocity
            dFramedt = gradient(double(trajData_r.frame_number))';
            midbodySpeed = sqrt(dmidbody_xdt.^2 + dmidbody_ydt.^2)./dFramedt;
            midbodyVelocity = [dmidbody_xdt; dmidbody_ydt]./dFramedt;
            % direction of segments pointing along midbody
            [~, dmidbody_yds] = gradient(squeeze(skelData_r(2,midbodyIndcs,:)),-1);
            [~, dmidbody_xds] = gradient(squeeze(skelData_r(1,midbodyIndcs,:)),-1);
            % sign speed based on relative orientation of velocity to midbody
            midbodySpeedSigned = getSignedSpeed(midbodyVelocity,[mean(dmidbody_xds); mean(dmidbody_yds)]);
            % ignore first and last frames of each worm's track
            wormChangeIndcs = gradient(double(trajData_r.worm_index_joined))~=0;
            midbodySpeedSigned(wormChangeIndcs)=NaN;
            % ignore frames with bad skeletonization
            midbodySpeedSigned(trajData_r.is_good_skel~=1)=NaN;
            % ignore skeletons otherwise filtered out
            midbodySpeedSigned(~trajData_r.filtered) = NaN;
            % smooth speed to denoise
            midbodySpeedSigned = smooth(midbodySpeedSigned,3,'moving');
            % find reversals in midbody speed
            [revStartInd, revDuration] = findReversals(...
                midbodySpeedSigned,trajData_r.worm_index_joined);
            % detect cluster status of reversals and features
            loneReversals = ismember(revStartInd,find(loneWorms));
            inClusterReversals = ismember(revStartInd,find(inCluster));
            smallClusterReversals = ismember(revStartInd,find(smallCluster));
            % estimate reversal frequencies and durations
            interRevTimes = diff(revStartInd);
            interRevTimesLone = diff(revStartInd(loneReversals));
            interRevTimesInCluster = diff(revStartInd(inClusterReversals));
            interRevTimesSmallCluster = diff(revStartInd(smallClusterReversals));
            % ignore those inter-reversal times that arise from
            % non-contiguous reversal sequences
            interRevTimesLone(diff(find(loneReversals))~=1) = NaN;
            interRevTimesInCluster(diff(find(inClusterReversals))~=1) = NaN;
            interRevTimesSmallCluster(diff(find(smallClusterReversals))~=1) = NaN;
            revDurationLone = revDuration(loneReversals);
            revDurationInCluster = revDuration(inClusterReversals);
            revDurationSmallCluster = revDuration(smallClusterReversals);
            % subtracting revDuration will more accurately reflect the
            % inter reversal time, and also set any data to NaN where
            % worm_index changed (where revDuration == NaN)
            interrevT_lone{fileCtr} = (interRevTimesLone - revDurationLone(1:end-1))/frameRate;
            interrevT_inCluster{fileCtr} = (interRevTimesInCluster - revDurationInCluster(1:end-1))/frameRate;
            interrevT_smallCluster{fileCtr} = (interRevTimesSmallCluster - revDurationSmallCluster(1:end-1))/frameRate;
            % counting reversal events
            Nrev_lone = nnz(loneReversals);
            Nrev_inCluster = nnz(inClusterReversals);
            Nrev_smallCluster = nnz(smallClusterReversals);
            T_lone = nnz(loneWorms)/frameRate;
            T_inCluster = nnz(inCluster)/frameRate;
            T_smallCluster = nnz(smallCluster)/frameRate;
            Trev_lone = nnz(midbodySpeedSigned(loneWorms)<0)/frameRate;
            Trev_inCluster = nnz(midbodySpeedSigned(inCluster)<0)/frameRate;
            Trev_smallCluster = nnz(midbodySpeedSigned(smallCluster)<0)/frameRate;
            reversalfreq_lone(fileCtr) = Nrev_lone./(T_lone - Trev_lone);
            reversalfreq_inCluster(fileCtr) = Nrev_inCluster./(T_inCluster - Trev_inCluster);
            reversalfreq_smallCluster(fileCtr) = Nrev_smallCluster./(T_smallCluster - Trev_smallCluster);
            reversaldurations_lone{fileCtr} = revDuration(loneReversals)/frameRate;
            reversaldurations_inCluster{fileCtr} = revDuration(inClusterReversals)/frameRate;
            reversaldurations_smallCluster{fileCtr} = revDuration(smallClusterReversals)/frameRate;
        end
        %pool data from all files
        interrevT_inCluster = vertcat(interrevT_inCluster{:});
        interrevT_lone = vertcat(interrevT_lone{:});
        interrevT_smallCluster = vertcat(interrevT_smallCluster{:});
        %% plot data
        % inter-reversal time
        set(0,'CurrentFigure',revInterTimeFig)
        ecdf(interrevT_lone,'Bounds','on','function','survivor')
        hold on
        if ~strcmp(wormnum,'1W')
            if length(interrevT_smallCluster)~=1
                ecdf(interrevT_smallCluster,'Bounds','on','function','survivor')
            end
            ecdf(interrevT_inCluster,'Bounds','on','function','survivor')
        end
        set(revInterTimeFig.Children,'YScale','log')
        title(revInterTimeFig.Children,[strains{strainCtr} ' ' wormnum],'FontWeight','normal');
        set(revInterTimeFig,'PaperUnits','centimeters')
        revInterTimeFig.Children.XLabel.String = 'inter-reversal time (s)';
        revInterTimeFig.Children.YLabel.String = 'cumulative probability';
        revInterTimeFig.Children.YLim(1) = 1e-2;
        revInterTimeFig.Children.XLim(2) = 120;
        if ~strcmp(wormnum,'1W')
            if length(interrevT_smallCluster)~=1
                legend(revInterTimeFig.Children.Children([9 6 3]),{'lone worms','small cluster','in cluster',})
            else
                legend(revInterTimeFig.Children.Children([6 3]),{'lone worms','in cluster'})
            end
        else
            legend(revInterTimeFig.Children,'single worms')
        end
        figurename = ['figures/reversals/reversalintertime_bodywall_' strains{strainCtr} '_' wormnum];
        exportfig(revInterTimeFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        
        % reversal frequency
        set(0,'CurrentFigure',revFreqFig)
        notBoxPlot([reversalfreq_lone,reversalfreq_smallCluster,reversalfreq_inCluster],...
            numCtr+[-0.3 0 0.3],'markMedian',true,'jitter',0.2)%,'style','line')

        revFreqFig.Children.XLim = [0 length(wormnums)+1];
        
        % reversal duration
        reversaldurations_lone = vertcat(reversaldurations_lone{:});
        reversaldurations_smallCluster = vertcat(reversaldurations_smallCluster{:});
        reversaldurations_inCluster = vertcat(reversaldurations_inCluster{:});
        histogram(revDurFig.Children,reversaldurations_lone,0:1/frameRate:15,...
            'Normalization','pdf','DisplayStyle','stairs');
        histogram(revDurFig.Children,reversaldurations_smallCluster,0:1/frameRate:15,...
            'Normalization','pdf','DisplayStyle','stairs','EdgeColor',0.5*ones(1,3));
        histogram(revDurFig.Children,reversaldurations_inCluster,0:1/frameRate:15,...
            'Normalization','pdf','DisplayStyle','stairs','EdgeColor','r');
        title(revDurFig.Children,[strains{strainCtr} ' ' wormnum],'FontWeight','normal');
        set(revDurFig,'PaperUnits','centimeters')
        xlabel(revDurFig.Children,'time (s)')
        ylabel(revDurFig.Children,'P')
        revDurFig.Children.XLim = [0 15];
        if ~strcmp(wormnum,'1W')
            legend(revDurFig.Children,{'lone worms','small cluster','in cluster'})
        else
            legend(revDurFig.Children,'single worms')
        end
        figurename = ['figures/reversals/reversaldurations_bodywall_' strains{strainCtr} '_' wormnum];
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
    figurename = ['figures/reversals/reversalfrequency_bodywall_' strains{strainCtr}];
    exportfig(revFreqFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end