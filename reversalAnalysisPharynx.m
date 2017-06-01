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
        for fileCtr = 1:numFiles % can be parfor?
            filename_g = filenames_g{fileCtr};
            trajData_g = h5read(filename_g,'/trajectories_data');
            blobFeats_g = h5read(filename_g,'/blob_features');
            skelData_g = h5read(filename_g,'/skeleton');
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
                loneWorms = min_neighbr_dist>=minNeighbrDist;
                inCluster = num_close_neighbrs>=3;
                smallCluster = (num_close_neighbrs==2 & neighbr_dist(:,3)>=minNeighbrDist)...
                    |(num_close_neighbrs==3 & neighbr_dist(:,4)>=minNeighbrDist)...
                    |(num_close_neighbrs==4 & neighbr_dist(:,5)>=minNeighbrDist);
            else
                loneWorms = true(size(trajData_g.frame_number));
                inCluster = false(size(trajData_g.frame_number));
                smallCluster = false(size(trajData_g.frame_number));
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
            [revStartInd, revDuration] = findReversals(...
                speedSigned,trajData_g.worm_index_joined);
            % detect cluster status of reversals and features
            loneReversals = ismember(revStartInd,find(loneWorms));
            inclusterReversals = ismember(revStartInd,find(inCluster));
            smallClusterReversals = ismember(revStartInd,find(smallCluster));
            % estimate reversal frequencies and durations
            interRevTimes = diff(revStartInd);
            interRevTimesLone = diff(revStartInd(loneReversals));
            interRevTimesCluster = diff(revStartInd(inclusterReversals));
            interRevTimesSmallCluster = diff(revStartInd(smallClusterReversals));
            % ignore those inter-reversal times that arise from
            % non-contiguous reversal sequences
            interRevTimesLone(diff(find(loneReversals))~=1) = NaN;
            interRevTimesCluster(diff(find(inclusterReversals))~=1) = NaN;
            interRevTimesSmallCluster(diff(find(smallClusterReversals))~=1) = NaN;
            revDurationLone = revDuration(loneReversals);
            revDurationCluster = revDuration(inclusterReversals);
            revDurationSmallCluster = revDuration(smallClusterReversals);
            % subtracting revDuration will more accurately reflect the
            % inter reversal time, and also set any data to NaN where
            % worm_index changed (where revDuration == NaN)
            interrevT_lone{fileCtr} = (interRevTimesLone - revDurationLone(1:end-1))/frameRate;
            interrevT_incluster{fileCtr} = (interRevTimesCluster - revDurationCluster(1:end-1))/frameRate;
            interrevT_smallCluster{fileCtr} = (interRevTimesSmallCluster - revDurationSmallCluster(1:end-1))/frameRate;
            % counting reversal events
            Nrev_lone = nnz(loneReversals);
            Nrev_incluster = nnz(inclusterReversals);
            Nrev_smallCluster = nnz(smallClusterReversals);
            T_lone = nnz(loneWorms)/frameRate;
            T_incluster = nnz(inCluster)/frameRate;
            T_smallcluster = nnz(smallCluster)/frameRate;
            Trev_lone = nnz(speedSigned(loneWorms)<0)/frameRate;
            Trev_incluster = nnz(speedSigned(inCluster)<0)/frameRate;
            Trev_smallcluster = nnz(speedSigned(smallClusterReversals)<0)/frameRate;
            reversalfreq_lone(fileCtr) = Nrev_lone./(T_lone - Trev_lone);
            reversalfreq_incluster(fileCtr) = Nrev_incluster./(T_incluster - Trev_incluster);
            reversalfreq_smallcluster(fileCtr) = Nrev_smallCluster./(T_smallcluster - Trev_smallcluster);
            reversaldurations_lone{fileCtr} = revDuration(loneReversals)/frameRate;
            reversaldurations_incluster{fileCtr} = revDuration(inclusterReversals)/frameRate;
            reversaldurations_smallcluster{fileCtr} = revDuration(smallClusterReversals)/frameRate;
        end
        %pool data from all files
        interrevT_lone = vertcat(interrevT_lone{:});
        interrevT_incluster = vertcat(interrevT_incluster{:});
        interrevT_smallCluster = vertcat(interrevT_smallCluster{:});
        %% plot data
        set(0,'CurrentFigure',revInterTimeFig)
        ecdf(interrevT_lone,'Bounds','on','function','survivor')
        hold on
        if ~strcmp(wormnum,'1W')
            ecdf(interrevT_smallCluster,'Bounds','on','function','survivor')
            ecdf(interrevT_incluster,'Bounds','on','function','survivor')
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
        histogram(revDurFig.Children,reversaldurations_lone,0:1/frameRate:15,...
            'Normalization','pdf','DisplayStyle','stairs');
        histogram(revDurFig.Children,reversaldurations_smallcluster,0:1/frameRate:15,...
            'Normalization','pdf','DisplayStyle','stairs','EdgeColor',0.5*ones(1,3));
        histogram(revDurFig.Children,reversaldurations_incluster,0:1/frameRate:15,...
            'Normalization','pdf','DisplayStyle','stairs','EdgeColor','r');
        %
        title(revDurFig.Children,[strains{strainCtr} ' ' wormnum],'FontWeight','normal');
        set(revDurFig,'PaperUnits','centimeters')
        xlabel(revDurFig.Children,'time (s)')
        ylabel(revDurFig.Children,'P')
        revDurFig.Children.XLim = [0 10];
        if ~strcmp(wormnum,'1W')
            legend(revDurFig.Children,{'lone worms','small cluster','in cluster'})
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