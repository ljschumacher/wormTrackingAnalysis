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

% set parameters for analysis
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
dataset = 1; % 1 or 2

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
if dataset ==1
    strains = {'N2','HA','npr1'};
elseif dataset ==2
    strains = {'N2','npr1'};
end

% analysis
for strainCtr = 1:length(strains)
    revFreqFig = figure; hold on
    for numCtr = 1:length(wormnums)
        revDurFig = figure; hold on
        revInterTimeFig = figure; hold on
        wormnum = wormnums{numCtr};
        %% load data
        if dataset ==1
            filenames_g = importdata(['datalists/' strains{strainCtr} '_' wormnum '_list.txt']);
        elseif dataset ==2
            filenames_g = importdata(['datalists/' strains{strainCtr} '_' wormnum '_g_list.txt']);
        end
        numFiles = length(filenames_g);
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
            filename_g = filenames_g{fileCtr};
            trajData_g = h5read(filename_g,'/trajectories_data');
            blobFeats_g = h5read(filename_g,'/blob_features');
            skelData_g = h5read(filename_g,'/skeleton');
            frameRate = double(h5readatt(filename_g,'/plate_worms','expected_fps'));
            % filter green data by blob size and intensity
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
                smallCluster = trajData_g.filtered& trajData_g.has_skeleton & ...
                    ((num_close_neighbrs == 2 & neighbr_dist(:,3)>=(loneClusterRadius))...
                    |(num_close_neighbrs == 3 & neighbr_dist(:,4)>=(loneClusterRadius))...
                    |(num_close_neighbrs == 4 & neighbr_dist(:,5)>=(loneClusterRadius)));
            else
                loneWorms = true(size(trajData_g.frame_number));
                inCluster = false(size(trajData_g.frame_number));
                smallCluster = false(size(trajData_g.frame_number));
            end
            %% load signed speed from blobFeats
            % sign speed based on relative orientation of velocity to midbody
            speedSigned = blobFeats_g.signed_speed;
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
            Trev_lone = nnz(speedSigned(loneWorms)<0)/frameRate;
            Trev_inCluster = nnz(speedSigned(inCluster)<0)/frameRate;
            Trev_smallCluster = nnz(speedSigned(smallClusterReversals)<0)/frameRate;
            reversalfreq_lone(fileCtr) = Nrev_lone./(T_lone - Trev_lone);
            reversalfreq_inCluster(fileCtr) = Nrev_inCluster./(T_inCluster - Trev_inCluster);
            reversalfreq_smallCluster(fileCtr) = Nrev_smallCluster./(T_smallCluster - Trev_smallCluster);
            reversaldurations_lone{fileCtr} = revDuration(loneReversals)/frameRate;
            reversaldurations_inCluster{fileCtr} = revDuration(inClusterReversals)/frameRate;
            reversaldurations_smallCluster{fileCtr} = revDuration(smallClusterReversals)/frameRate;
        end
        %pool data from all files
        interrevT_lone = vertcat(interrevT_lone{:});
        interrevT_inCluster = vertcat(interrevT_inCluster{:});
        interrevT_smallCluster = vertcat(interrevT_smallCluster{:});
        interrevT_smallCluster(isnan(interrevT_smallCluster))=[];
        
        %% plot data
        % inter-reversal time
        set(0,'CurrentFigure',revInterTimeFig)
        ecdf(interrevT_lone,'Bounds','on','function','survivor')
        hold on
        if ~strcmp(wormnum,'1W')
            if length(interrevT_smallCluster)>=2
                ecdf(interrevT_smallCluster,'Bounds','on','function','survivor')
            end
            ecdf(interrevT_inCluster,'Bounds','on','function','survivor')
        end
        set(revInterTimeFig.Children,'YScale','log')
        title(revInterTimeFig.Children,[strains{strainCtr} ' ' wormnum],'FontWeight','normal');
        set(revInterTimeFig,'PaperUnits','centimeters')
        revInterTimeFig.Children.XLabel.String = 'inter-reversal time (s)';
        revInterTimeFig.Children.YLabel.String = 'cumulative probability';
        %revInterTimeFig.Children.YLim(1) = 1e-2;
        revInterTimeFig.Children.XLim(2) = 120;
        if ~strcmp(wormnum,'1W')
            if length(interrevT_smallCluster)>=2
                legend(revInterTimeFig.Children.Children([9 6 3]),{'lone worms','small cluster','in cluster',})
            else
                legend(revInterTimeFig.Children.Children([6 3]),{'lone worms','in cluster'})
            end
        else
            legend(revInterTimeFig.Children,'single worms')
        end
        if dataset ==1
            figurename = ['figures/reversals/reversalintertime_pharynx1_' strains{strainCtr} '_' wormnum];
        elseif dataset ==2
            figurename = ['figures/reversals/reversalintertime_pharynx_' strains{strainCtr} '_' wormnum];
        end
        savefig([figurename '.fig'])
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
        reversaldurations_inCluster = vertcat(reversaldurations_inCluster{:});
        reversaldurations_smallCluster = vertcat(reversaldurations_smallCluster{:});
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
        revDurFig.Children.XLim = [0 10];
        if ~strcmp(wormnum,'1W')
            legend(revDurFig.Children,{'lone worms','small cluster','in cluster'})
        else
            legend(revDurFig.Children,'single worms')
        end
        if dataset ==1
            figurename = ['figures/reversals/reversaldurations_pharynx1_' strains{strainCtr} '_' wormnum];
        elseif dataset ==2
            figurename = ['figures/reversals/reversaldurations_pharynx_' strains{strainCtr} '_' wormnum];
        end
        savefig([figurename '.fig'])
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
    if dataset == 1
        figurename = ['figures/reversals/reversalfrequency_pharynx1_' strains{strainCtr}];
    elseif dataset ==2
        figurename = ['figures/reversals/reversalfrequency_pharynx_' strains{strainCtr}];
    end
    savefig([figurename '.fig'])
    exportfig(revFreqFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end