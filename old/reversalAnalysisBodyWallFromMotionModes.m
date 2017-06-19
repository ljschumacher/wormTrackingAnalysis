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
minNeighbrDist = 1500;

plotColors = lines(length(wormnums));

for strainCtr = 1:length(strains)
    revFreqFig = figure; hold on
    for numCtr = 1:length(wormnums)
        revDurFig = figure; hold on
        wormnum = wormnums{numCtr};
        %% load data
        filenames_r = importdata(['datalists/' strains{strainCtr} '_' wormnum '_r_list.txt']);
        if ~strcmp(wormnum,'1W')
            filenames_g = importdata(['datalists/' strains{strainCtr} '_' wormnum '_g_list.txt']);
        else
            filenames_g = {};
        end
        numFiles = length(filenames_r);
        reversalfreq_lone = NaN(numFiles,1);
        reversaldurations_lone = cell(numFiles,1);
        reversalfreq_incluster = NaN(numFiles,1);
        reversaldurations_incluster = cell(numFiles,1);
        reversalfreq_neither = NaN(numFiles,1);
        reversaldurations_neither = cell(numFiles,1);
        for fileCtr = 1:numFiles % can be parfor?
            filename_r = filenames_r{fileCtr};
            trajData_r = h5read(filename_r,'/trajectories_data');
            blobFeats_r = h5read(filename_r,'/blob_features');
            skelData_r = h5read(filename_r,'/skeleton');
            if ~strcmp(wormnum,'1W')
                filename_g = filenames_g{fileCtr};
                trajData_g = h5read(filename_g,'/trajectories_data');
                % filter data
                blobFeats_g = h5read(filename_g,'/blob_features');
                trajData_g.filtered = filterIntensityAndSize(blobFeats_g,pixelsize,...
                    intensityThresholds_g(wormnum),maxBlobSize_g);
            end
            featData_r = h5read(strrep(filename_r,'skeletons','features'),'/features_timeseries');
            frameRate = double(h5readatt(filename_r,'/plate_worms','expected_fps'));
            % filter by blob size and intensity
            if contains(filename_r,'55')||contains(filename_r,'54')
                intensityThreshold_r = 80;
            else
                intensityThreshold_r = 40;
            end
            trajData_r.filtered = filterIntensityAndSize(blobFeats_r,pixelsize,...
                intensityThreshold_r,maxBlobSize);
            % filter by skeleton length
            trajData_r.filtered = trajData_r.filtered&logical(trajData_r.is_good_skel)&...
                filterSkelLength(skelData_r,pixelsize,minSkelLength,maxSkelLength);
            %% calculate stats
            if ~strcmp(wormnum,'1W')
                min_neighbr_dist = h5read(filename_r,'/min_neighbr_dist');
                num_close_neighbrs = h5read(filename_r,'/num_close_neighbrs');
                loneWorms = min_neighbr_dist>=minNeighbrDist;
                inCluster = num_close_neighbrs>=3;
                neitherClusterNorLone = num_close_neighbrs==1|num_close_neighbrs==2;
            else
                loneWorms = true(size(trajData_r.frame_number));
                inCluster = false(size(trajData_r.frame_number));
                neitherClusterNorLone = false(size(trajData_r.frame_number));
            end
            % features from the tracker
            % ignore data otherwise filtered out
            featData_r.motion_modes(ismember(featData_r.skeleton_id+1,find(~trajData_r.filtered))) = NaN;
            %%
            % find reversals in motion modes
            [revStartInd, revDuration] = findReversals(...
                featData_r.motion_modes,featData_r.worm_index);
            % detect cluster status of reversals and features
            loneReversals = ismember(featData_r.skeleton_id(revStartInd)+1,find(loneWorms));
            inclusterReversals = ismember(featData_r.skeleton_id(revStartInd)+1,find(inCluster));
            neitherClusterNorLoneReversals = ismember(featData_r.skeleton_id(revStartInd)+1,find(neitherClusterNorLone));
            loneWormsFeats = ismember(featData_r.skeleton_id+1,find(loneWorms));
            inClusterFeats = ismember(featData_r.skeleton_id+1,find(inCluster));
            neitherClusterNorLoneFeats = ismember(featData_r.skeleton_id+1,find(neitherClusterNorLone));
            % estimate reversal frequencies and durations
            Nrev_lone = nnz(loneReversals);
            Nrev_incluster = nnz(inclusterReversals);
            Nrev_neither = nnz(neitherClusterNorLoneReversals);
            T_lone = nnz(loneWormsFeats)/frameRate;
            T_incluster = nnz(inClusterFeats)/frameRate;
            T_neither = nnz(neitherClusterNorLone)/frameRate;
            Trev_lone = nnz(featData_r.motion_modes(loneWormsFeats)<0)/frameRate;
            Trev_incluster = nnz(featData_r.motion_modes(inClusterFeats)<0)/frameRate;
            Trev_neither = nnz(featData_r.motion_modes(neitherClusterNorLoneFeats)<0)/frameRate;
            reversalfreq_lone(fileCtr) = Nrev_lone./(T_lone - Trev_lone);
            reversalfreq_incluster(fileCtr) = Nrev_incluster./(T_incluster - Trev_incluster);
            reversalfreq_neither(fileCtr) = Nrev_neither./(T_neither - Trev_neither);
            reversaldurations_lone{fileCtr} = revDuration(loneReversals)/frameRate;
            reversaldurations_incluster{fileCtr} = revDuration(inclusterReversals)/frameRate;
            reversaldurations_neither{fileCtr} = revDuration(neitherClusterNorLoneReversals)/frameRate;
        end
        %% plot data
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
        figurename = ['figures/reversaldurations_bodywall_FromMotionModes_' strains{strainCtr} '_' wormnum];
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
    figurename = ['figures/reversalfrequency_bodywall_FromMotionModes_' strains{strainCtr}];
    exportfig(revFreqFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end