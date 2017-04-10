% calculate various single worm statistics for different numbers of worms
% on plate

% issues / todo:
% - BUILD IN MINIMUM SPEED
% THRESHOLD TO AVOID REVERSAL ARTEFACTS FROM NOISY TRACKING?
% - is it better to plot a ecdf for all data or a mean ecdf?
% - filter less for detection of reversals (include bad skels?)?
% - assign untracked reversal starts/ends by estimating from last +/- and
% first -/+ velocity?
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
intensityThresholds_g = [100, 60, 40];
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
                trajData_g.filtered = (blobFeats_g.area*pixelsize^2<=maxBlobSize_g)&...
                    (blobFeats_g.intensity_mean>=intensityThresholds_g(numCtr));
            end
            featData = h5read(strrep(filename,'skeletons','features'),'/features_timeseries');
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            % filter by blob size and intensity
            if contains(filename,'55')||contains(filename,'54')
                intensityThreshold = 70;
            else
                intensityThreshold = 35;
            end
            trajData.filtered = (blobFeats.area*pixelsize^2<=maxBlobSize)&...
                (blobFeats.intensity_mean>=intensityThreshold)&...
                logical(trajData.is_good_skel);
            % filter by skeleton length
            skelLengths = sum(sqrt(sum((diff(skelData,1,2)*pixelsize).^2)));
            trajData.filtered(skelLengths>maxSkelLength|skelLengths<minSkelLength)...
                = false;
            %% calculate stats
            if ~strcmp(wormnum,'1W')
                try
                    min_neighbor_dist_rr = h5read(filename,'/min_neighbor_dist_rr');
                    min_neighbor_dist_rg = h5read(filename,'/min_neighbor_dist_rg');
                    num_close_neighbours_rg = h5read(filename,'/num_close_neighbours_rg');
                catch
                    disp(['Calculating cluster status from ' filename ])
                    [min_neighbor_dist_rr, min_neighbor_dist_rg, num_close_neighbours_rg] ...
                        = calculateClusterStatus(trajData,trajData_g,pixelsize,500);
                    % write stats to hdf5-file
                    h5create(filename,'/min_neighbor_dist_rr',...
                        size(min_neighbor_dist_rr))
                    h5write(filename,'/min_neighbor_dist_rr',...
                        single(min_neighbor_dist_rr))
                    h5create(filename,'/min_neighbor_dist_rg',...
                        size(min_neighbor_dist_rg))
                    h5write(filename,'/min_neighbor_dist_rg',...
                        single(min_neighbor_dist_rg))
                    h5create(filename,'/num_close_neighbours_rg',...
                        size(num_close_neighbours_rg))
                    h5write(filename,'/num_close_neighbours_rg',...
                        uint16(num_close_neighbours_rg))
                end
                loneWorms = min_neighbor_dist_rr>=1100&min_neighbor_dist_rg>=1600;
                inCluster = num_close_neighbours_rg>=3;
                neitherClusterNorLone = num_close_neighbours_rg==1|num_close_neighbours_rg==2;
                %~inCluster&~loneWorms;
            else
                loneWorms = true(size(trajData.frame_number));
                inCluster = false(size(trajData.frame_number));
                neitherClusterNorLone = false(size(trajData.frame_number));
            end
            % features from the tracker
            midbodySpeedSigned = featData.midbody_speed;
            % ignore data otherwise filtered out
            midbodySpeedSigned(ismember(featData.skeleton_id+1,find(~trajData.filtered))) = NaN;
            % smooth speed to denoise
            midbodySpeedSigned = smooth(midbodySpeedSigned,3,'moving');
            %%
            % find reversals in midbody speed
            [revStartInd, revDuration] = findReversals(...
                midbodySpeedSigned,featData.worm_index);
            % detect cluster status of reversals and features
            loneReversals = ismember(featData.skeleton_id(revStartInd)+1,find(loneWorms));
            inclusterReversals = ismember(featData.skeleton_id(revStartInd)+1,find(inCluster));
            neitherClusterNorLoneReversals = ismember(featData.skeleton_id(revStartInd)+1,find(neitherClusterNorLone));
            loneWormsFeats = ismember(featData.skeleton_id+1,find(loneWorms));
            inClusterFeats = ismember(featData.skeleton_id+1,find(inCluster));
            neitherClusterNorLoneFeats = ismember(featData.skeleton_id+1,find(neitherClusterNorLone));
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
            T_lone = nnz(loneWormsFeats)/frameRate;
            T_incluster = nnz(inClusterFeats)/frameRate;
            T_neither = nnz(neitherClusterNorLone)/frameRate;
            Trev_lone = nnz(midbodySpeedSigned(loneWormsFeats)<0)/frameRate;
            Trev_incluster = nnz(midbodySpeedSigned(inClusterFeats)<0)/frameRate;
            Trev_neither = nnz(midbodySpeedSigned(neitherClusterNorLoneFeats)<0)/frameRate;
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
            legend(revInterTimeFig.Children.Children([9 6 3]),{'lone worms','in cluster','neither'})
        else
            legend(revInterTimeFig.Children,'single worms')
        end
        figurename = ['figures/reversalintertimeFromFeatures_' strains{strainCtr} '_' wormnum];
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
        figurename = ['figures/reversaldurationsFromFeatures_' strains{strainCtr} '_' wormnum];
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
    figurename = ['figures/reversalfrequencyFromFeatures_' strains{strainCtr}];
    exportfig(revFreqFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end