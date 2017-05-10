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
maxBlobSize = 2.5e5;
maxBlobSize_g = 1e4;
minSkelLength = 850;
maxSkelLength = 1500;
midbodyIndcs = 19:33;
plotColors = lines(length(wormnums));

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
        reversalfreq_incluster = NaN(numFiles,1);
        reversaldurations_incluster = cell(numFiles,1);
        reversalfreq_neither = NaN(numFiles,1);
        reversaldurations_neither = cell(numFiles,1);
        interrevT_lone = cell(numFiles,1);
        interrevT_incluster = cell(numFiles,1);
        interrevT_neither = cell(numFiles,1);
        for fileCtr = 1:numFiles % can be parfor?
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton');
            if ~strcmp(wormnum,'1W')
                filename_g = filenames_g{fileCtr};
                trajData_g = h5read(filename_g,'/trajectories_data');
                % filter data
                blobFeats_g = h5read(filename_g,'/blob_features');
                trajData_g.filtered = filterIntensityAndSize(blobFeats_g,pixelsize,...
                    intensityThresholds_g(wormnum),maxBlobSize_g);
            end
            % %             featData = h5read(strrep(filename,'skeletons','features'),'/features_timeseries');
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            % filter by blob size and intensity
            if contains(filename,'55')||contains(filename,'54')
                intensityThreshold_r = 80;
            else
                intensityThreshold_r = 40;
            end
            trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                    intensityThreshold_r,maxBlobSize);
            % filter by skeleton length
            trajData.filtered = trajData.filtered&logical(trajData.is_good_skel)&...
                filterSkelLength(skelData,pixelsize,minSkelLength,maxSkelLength);
            %% calculate stats
            if ~strcmp(wormnum,'1W')
                min_neighbr_dist_rr = h5read(filename,'/min_neighbr_dist_rr');
                min_neighbr_dist_rg = h5read(filename,'/min_neighbr_dist_rg');
                num_close_neighbrs_rg = h5read(filename,'/num_close_neighbrs_rg');
                loneWorms = min_neighbr_dist_rr>=1100&min_neighbr_dist_rg>=1600;
                inCluster = num_close_neighbrs_rg>=3;
                neitherClusterNorLone = num_close_neighbrs_rg==1|num_close_neighbrs_rg==2;
            else
                loneWorms = true(size(trajData.frame_number));
                inCluster = false(size(trajData.frame_number));
                neitherClusterNorLone = false(size(trajData.frame_number));
            end
            %% calculate overall midbody speeds, until we can link trajectories and features from the tracker
            % centroids of midbody skeleton
            midbody_x = mean(squeeze(skelData(1,midbodyIndcs,:)))*pixelsize;
            midbody_y = mean(squeeze(skelData(2,midbodyIndcs,:)))*pixelsize;
            % change in centroid position over time (issue of worm shifts?)
            dmidbody_xdt = gradient(midbody_x)*frameRate;
            dmidbody_ydt = gradient(midbody_y)*frameRate;
            % midbody speed and velocity
            midbodySpeed = sqrt(dmidbody_xdt.^2 + dmidbody_ydt.^2)./gradient(double(trajData.frame_number))';
            midbodyVelocity = [dmidbody_xdt; dmidbody_ydt]./gradient(double(trajData.frame_number))';
            % direction of segments pointing along midbody
            [~, dmidbody_yds] = gradient(squeeze(skelData(2,midbodyIndcs,:)),-1);
            [~, dmidbody_xds] = gradient(squeeze(skelData(1,midbodyIndcs,:)),-1);
            % sign speed based on relative orientation of velocity to midbody
            midbodySpeedSigned = sign(sum(midbodyVelocity.*[mean(dmidbody_xds); mean(dmidbody_yds)])).*midbodySpeed;
            % ignore first and last frames of each worm's track
            wormChangeIndcs = gradient(double(trajData.worm_index_joined))~=0;
            midbodySpeedSigned(wormChangeIndcs)=NaN;
            % ignore frames with bad skeletonization
            midbodySpeedSigned(trajData.is_good_skel~=1)=NaN;
            % ignore skeletons otherwise filtered out
            midbodySpeedSigned(~trajData.filtered) = NaN;
            % smooth speed to denoise
            midbodySpeedSigned = smooth(midbodySpeedSigned,3,'moving');
            % find reversals in midbody speed
            [revStartInd, revDuration] = findReversals(...
                midbodySpeedSigned,trajData.worm_index_joined);
            % detect cluster status of reversals and features
            loneReversals = ismember(revStartInd,find(loneWorms));
            inclusterReversals = ismember(revStartInd,find(inCluster));
            neitherClusterNorLoneReversals = ismember(revStartInd,find(neitherClusterNorLone));
            % estimate reversal frequencies and durations
            interRevTimes = diff(revStartInd);
            interRevTimesLone = diff(revStartInd(loneReversals));
            interRevTimesCluster = diff(revStartInd(inclusterReversals));
            interRevTimesNeither = diff(revStartInd(neitherClusterNorLoneReversals));
            % ignore those inter-reversal times that arise from
            % non-contiguous reversal sequences
            interRevTimesLone(diff(find(loneReversals))~=1) = NaN;
            interRevTimesCluster(diff(find(inclusterReversals))~=1) = NaN;
            interRevTimesNeither(diff(find(neitherClusterNorLoneReversals))~=1) = NaN;
            revDurationLone = revDuration(loneReversals);
            revDurationCluster = revDuration(inclusterReversals);
            revDurationNeither = revDuration(neitherClusterNorLoneReversals);
            % subtracting revDuration will more accurately reflect the
            % inter reversal time, and also set any data to NaN where
            % worm_index changed (where revDuration == NaN)
            interrevT_lone{fileCtr} = (interRevTimesLone - revDurationLone(1:end-1))/frameRate;
            interrevT_incluster{fileCtr} = (interRevTimesCluster - revDurationCluster(1:end-1))/frameRate;
            interrevT_neither{fileCtr} = (interRevTimesNeither - revDurationNeither(1:end-1))/frameRate;
            % counting reversal events
            Nrev_lone = nnz(loneReversals);
            Nrev_incluster = nnz(inclusterReversals);
            Nrev_neither = nnz(neitherClusterNorLoneReversals);
            T_lone = nnz(loneWorms)/frameRate;
            T_incluster = nnz(inCluster)/frameRate;
            T_neither = nnz(neitherClusterNorLone)/frameRate;
            Trev_lone = nnz(midbodySpeedSigned(loneWorms)<0)/frameRate;
            Trev_incluster = nnz(midbodySpeedSigned(inCluster)<0)/frameRate;
            Trev_neither = nnz(midbodySpeedSigned(neitherClusterNorLoneReversals)<0)/frameRate;
            reversalfreq_lone(fileCtr) = Nrev_lone./(T_lone - Trev_lone);
            reversalfreq_incluster(fileCtr) = Nrev_incluster./(T_incluster - Trev_incluster);
            reversalfreq_neither(fileCtr) = Nrev_neither./(T_neither - Trev_neither);
            reversaldurations_lone{fileCtr} = revDuration(loneReversals)/frameRate;
            reversaldurations_incluster{fileCtr} = revDuration(inclusterReversals)/frameRate;
            reversaldurations_neither{fileCtr} = revDuration(neitherClusterNorLoneReversals)/frameRate;
        end
        %pool data from all files
        interrevT_lone = vertcat(interrevT_lone{:});
        interrevT_incluster = vertcat(interrevT_incluster{:});
        interrevT_neither = vertcat(interrevT_neither{:});
        %% plot data
        set(0,'CurrentFigure',revInterTimeFig)
        ecdf(interrevT_lone,'Bounds','on','function','survivor')
        hold on
        if ~strcmp(wormnum,'1W')
            ecdf(interrevT_incluster,'Bounds','on','function','survivor')
            ecdf(interrevT_neither,'Bounds','on','function','survivor')
        end
        set(revInterTimeFig.Children,'YScale','log')
        title(revInterTimeFig.Children,[strains{strainCtr} ' ' wormnum],'FontWeight','normal');
        set(revInterTimeFig,'PaperUnits','centimeters')
        revInterTimeFig.Children.XLabel.String = 'inter-reversal time (s)';
        revInterTimeFig.Children.YLabel.String = 'cumulative probability';
        revInterTimeFig.Children.YLim(1) = 1e-2;
        revInterTimeFig.Children.XLim(2) = 120;
        if ~strcmp(wormnum,'1W')
            legend(revInterTimeFig.Children.Children([9 6 3]),{'lone worms','in cluster','neither',})
        else
            legend(revInterTimeFig.Children,'single worms')
        end
        figurename = ['figures/reversalintertime_' strains{strainCtr} '_' wormnum];
        exportfig(revInterTimeFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        %
        set(0,'CurrentFigure',revFreqFig)
        notBoxPlot([reversalfreq_lone,reversalfreq_neither,reversalfreq_incluster],...
            numCtr+[-0.3 0 0.3],'markMedian',true,'jitter',0.2)%,'style','line')
        %         boxplot(revFreqFig.Children,reversalfreq_lone,'Positions',numCtr-1/4,...
        %             'Notch','off')
        %         boxplot(revFreqFig.Children,reversalfreq_neither,'Positions',numCtr,...
        %             'Notch','off','Colors',0.5*ones(1,3))
        %         boxplot(revFreqFig.Children,reversalfreq_incluster,'Positions',numCtr+1/4,...
        %             'Notch','off','Colors','r')
        revFreqFig.Children.XLim = [0 length(wormnums)+1];
        %
        reversaldurations_lone = vertcat(reversaldurations_lone{:});
        reversaldurations_incluster = vertcat(reversaldurations_incluster{:});
        reversaldurations_neither = vertcat(reversaldurations_neither{:});
        histogram(revDurFig.Children,reversaldurations_lone,0:1/frameRate:15,...
            'Normalization','pdf','DisplayStyle','stairs');
        histogram(revDurFig.Children,reversaldurations_neither,0:1/frameRate:15,...
            'Normalization','pdf','DisplayStyle','stairs','EdgeColor',0.5*ones(1,3));
        histogram(revDurFig.Children,reversaldurations_incluster,0:1/frameRate:15,...
            'Normalization','pdf','DisplayStyle','stairs','EdgeColor','r');
        %
        title(revDurFig.Children,[strains{strainCtr} ' ' wormnum],'FontWeight','normal');
        set(revDurFig,'PaperUnits','centimeters')
        xlabel(revDurFig.Children,'time (s)')
        ylabel(revDurFig.Children,'P')
        revDurFig.Children.XLim = [0 15];
        if ~strcmp(wormnum,'1W')
            legend(revDurFig.Children,{'lone worms','neither','in cluster'})
        else
            legend(revDurFig.Children,'single worms')
        end
        figurename = ['figures/reversaldurations_' strains{strainCtr} '_' wormnum];
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
    figurename = ['figures/reversalfrequency_' strains{strainCtr}];
    exportfig(revFreqFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end