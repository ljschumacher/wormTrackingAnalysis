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
load ~/Dropbox/Utilities/colormaps_ascii/increasing_cool/cmap_Blues.txt

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
useJoinedTraj = true
maxBlobSize_g = 1e4;
minPathLength = 50; % minimum path length of reversals to be included
knnbrNumValues =4:7; % which numbers of knnbrs to consider for checking density dependence
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
speedBinLimits = containers.Map({'npr1','N2'},{400, 300});
densityBinLimits = containers.Map({'npr1','N2'},{10, 5});
%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    %% load data
    if dataset ==1
        [phaseFrames,filenames,~] = xlsread(['datalists/' strain '_' wormnum '_list.xlsx'],1,'A1:E15','basic');
    elseif dataset ==2
        [phaseFrames,filenames,~] = xlsread(['datalists/' strain '_' wormnum '_g_list.xlsx'],1,'A1:E15','basic');
    end
    if ~useJoinedTraj
        filenames = strrep(filenames,'/data2/shared/data/twoColour/Results/',...
            '/end/home/lschumac/databackup/data/twoColour/ResultsUnjoinedTrajectories/');
    end
    %% intialize variables
    numFiles = length(filenames);
    knndensity = cell(numFiles,1);
    knndensityChange = cell(numFiles,1);
    knndensityAtRev = cell(numFiles,1);
    knndensityAtFwd = cell(numFiles,1);
    knndensityChangeAtRev = cell(numFiles,1);
    knndensityChangeAtFwd = cell(numFiles,1);
    signedSpeeds = cell(numFiles,1);
    %% loop through files and load data
    for fileCtr = 1:numFiles % can be parfor
        filename = filenames{fileCtr};
        shortFileName = filename(end-31:end-5);
        trajData = h5read(filename,'/trajectories_data');
        blobFeats = h5read(filename,'/blob_features');
        skelData = h5read(filename,'/skeleton');
        assert(size(skelData,1)==2&&size(skelData,2)==2,['unexpected skeleton size for ' filename]);
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        if frameRate == 0
            warning(['frame rate is zero for ' filename])
        end
        % filter by blob size and intensity
        trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
            intensityThresholds_g(wormnum),maxBlobSize_g);
        trajData.has_skeleton = squeeze(~any(any(isnan(skelData)))); % reset skeleton flag for pharynx data
        % check worm-indices are monotonically increasing
        assert(~any(diff(trajData.worm_index_joined)<0),['worm indices are not sorted as expected for ' filename])
        %% get relevant experimental phase
        if ~strcmp(wormnum,'1W')
            % apply phase restriction
            [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr);
            phaseFrameLogInd = trajData.frame_number < lastFrame & trajData.frame_number > firstFrame;
        else
            phaseFrameLogInd = true(size(trajData.frame_number));
        end
        %% calculate stats
        if any(phaseFrameLogInd)
            %% sign speed based on relative orientation of velocity to midbody
            signedSpeedThisFile = blobFeats.signed_speed*pixelsize*frameRate;
            numFramesTotal = length(signedSpeedThisFile);
            disp([shortFileName ' has ' num2str(100*nnz(isnan(signedSpeedThisFile))/numFramesTotal,2) '% of speeds = NaN'])
            % ignore first and last frames of each worm's track
            assert(all(trajData.worm_index_joined(2:end)>=trajData.worm_index_joined(1:end-1)),'worm indices are not sorted in monotonically increasing order')
            wormChangeIndcs = gradient(double(trajData.worm_index_joined))~=0;
            disp(['excluding ' num2str(100*nnz(wormChangeIndcs)/numFramesTotal,2) '% of speeds due to change in worm index'])
            signedSpeedThisFile(wormChangeIndcs)=NaN;
            % ignore frames without skeletonization
            disp(['excluding ' num2str(100*nnz(~trajData.has_skeleton)/numFramesTotal,2) '% of speeds due to missing skeletonization'])
            signedSpeedThisFile(~trajData.has_skeleton)=NaN;
            % ignore skeletons otherwise filtered out
            disp(['excluding ' num2str(100*nnz(~trajData.filtered)/numFramesTotal,2) '% of speeds due to other filtering'])
            signedSpeedThisFile(~trajData.filtered) = NaN;
            disp([shortFileName ' now has ' num2str(100*nnz(isnan(signedSpeedThisFile))/numFramesTotal) '% of speeds = NaN'])
            % smooth speed to denoise
            signedSpeedThisFile = smoothdata(signedSpeedThisFile,'movmean',round(frameRate/2),'omitnan');
            % ignore frames outside of specified phase
            if ~strcmp(phase,'fullMovie')
                signedSpeedThisFile(~phaseFrameLogInd) =NaN;
            end
            nonNanFramesLogInd = ~isnan(signedSpeedThisFile);
            signedSpeeds{fileCtr} = signedSpeedThisFile(nonNanFramesLogInd);
            %% find reversals in signed speed
            [revStartInd, revDuration, untrackedRevEnds, interRevTime, incompleteInterRev] = ...
                findReversals(signedSpeedThisFile,trajData.worm_index_joined,minPathLength,frameRate);
            [fwdStartInd, ~, ~, ~, ~] = ...
                findReversals(-signedSpeedThisFile,trajData.worm_index_joined,minPathLength,frameRate);
%             fwdStartInd = find(signedSpeedThisFile(1:end-1)<0&signedSpeedThisFile(2:end)>0);
            disp([num2str(100*mean(incompleteInterRev),2) '% of tracks are lost between reversals'])
            disp([num2str(100*mean(untrackedRevEnds),2) '% of tracks are lost during reversals'])
            % if we subtract rev duration from interrevtime (below), we
            % need to set all reversals with untracked ends as incomplete interrevs
            incompleteInterRev = incompleteInterRev|untrackedRevEnds;
            disp([num2str(100*mean(incompleteInterRev),2) '% of reversals are thus not fully tracked' newline])
            %% estimate density based on k-nearest neighbours
            neighbr_dist = h5read(filename,'/neighbr_distances')/1000;% convert to mm
            knndensityThisFile = zeros(numFramesTotal,10);
            for k = 1:10
                knndensityThisFile(:,k) = k./(pi*neighbr_dist(:,k).^2);
            end
            knndensity{fileCtr} = knndensityThisFile(nonNanFramesLogInd,:);            
            %% store reversal statistics for this file
            % densities at reversal
            knndensityAtRev{fileCtr} = knndensityThisFile(revStartInd,:);
            knndensityAtFwd{fileCtr} = knndensityThisFile(fwdStartInd,:);
                        %% estimate change in density
            % ignore frames where worm index changes before gradient calc
            knndensityThisFile(wormChangeIndcs) = NaN;
            % smooth density for 1s
            knndensityThisFile = smoothdata(knndensityThisFile,1,'movmean',round(frameRate),'omitnan');
            knndensityChangeThisFile = gradient(knndensityThisFile')'*frameRate; % change per second
            knndensityChange{fileCtr} = knndensityChangeThisFile(nonNanFramesLogInd,:);
            knndensityChangeAtRev{fileCtr} = knndensityChangeThisFile(revStartInd,:);
            knndensityChangeAtFwd{fileCtr} = knndensityChangeThisFile(fwdStartInd,:);
        end
    end
    %% pool data from all files
    knndensity = vertcat(knndensity{:});
    knndensityChange = vertcat(knndensityChange{:});
    knndensityAtRev = vertcat(knndensityAtRev{:});
    knndensityAtFwd = vertcat(knndensityAtFwd{:});
    knndensityChangeAtRev = vertcat(knndensityChangeAtRev{:});
    knndensityChangeAtFwd = vertcat(knndensityChangeAtFwd{:});
    signedSpeeds = vertcat(signedSpeeds{:});
    densityFig = figure; hold on
    speedAbsHist2Dfig = figure; hold on, hold off % initialize axes
    colormap(speedAbsHist2Dfig,flipud(cmap_Blues))
    speedSignedHist2Dfig = figure; hold on, hold off % initialize axes
    colormap(speedSignedHist2Dfig,flipud(cmap_Blues))
    revFig = figure; hold on, hold off % initialize axes
    fwdFig = figure; hold on, hold off % initialize axes
    revfwdFig = figure; hold on, hold off % initialize axes
    revFig2 = figure; hold on, hold off % initialize axes
    fwdFig2 = figure; hold on, hold off % initialize axes
    revfwdFig2 = figure; hold on, hold off % initialize axes
    for k = knnbrNumValues
        %% plot "convergence" of knn density estimate
        histogram(densityFig.Children, knndensity(:,k),'Normalization','Probability','DisplayStyle','stairs','BinLimits',[0 20])
        %% plot histogram of speeds vs density
        % for signed speed
        h = histogram2(speedSignedHist2Dfig.Children,knndensity(:,k),signedSpeeds,...
            'DisplayStyle','tile','YBinLimits',[-1 1].*speedBinLimits(strain),...
            'XBinLimits',[0 densityBinLimits(strain)],'EdgeColor','none','Normalization','Probability');
        normfactor = sum(h.BinCounts,2); % for conditional normalisation
        normfactor(normfactor==0) = 1;
        h.BinCounts = h.BinCounts./normfactor; % conditional normalisation 
        colormap(flipud(cmap_Blues))
        xlabel(speedSignedHist2Dfig.Children,['\rho_' num2str(k) ' (worms/mm^2)'])
        ylabel(speedSignedHist2Dfig.Children,'signed speed')
        figurename = ['figures/reversals/phaseSpecific/speedsigned-density_knn' num2str(k) ...
            '_pharynx_' strain '_' wormnum '_' phase '_data' num2str(dataset)];
        formatAndExportFigure(speedSignedHist2Dfig,figurename,useJoinedTraj,exportOptions)
        % for absolute speed
        h = histogram2(speedAbsHist2Dfig.Children,knndensity(:,k),abs(signedSpeeds),...
            'DisplayStyle','tile','YBinLimits',[0 speedBinLimits(strain)],...
            'XBinLimits',[0 densityBinLimits(strain)],'EdgeColor','none','Normalization','Probability');
        normfactor = sum(h.BinCounts,2); % for conditional normalisation
        normfactor(normfactor==0) = 1;
        h.BinCounts = h.BinCounts./normfactor; % conditional normalisation
        xlabel(speedAbsHist2Dfig.Children,['\rho_' num2str(k) ' (worms/mm^2)'])
        ylabel(speedAbsHist2Dfig.Children,'speed')
        figurename = ['figures/reversals/phaseSpecific/speedabs-density_knn' num2str(k) ...
            '_pharynx_' strain '_' wormnum '_' phase '_data' num2str(dataset)];
        formatAndExportFigure(speedAbsHist2Dfig,figurename,useJoinedTraj,exportOptions)
        %% plot histogram of increase in reversal frequency vs density
        [~, densitybins] = histcounts(knndensityAtRev(:,k),'BinLimits',[0 densityBinLimits(strain)]);
        % subsample to estimate variability
        numSubSamples = 100;
        sampleSizeFraction = 0.05;
        revRateSubsamples = subsampleRevFreq(knndensityAtRev(:,k),knndensity(:,k),...
            densitybins,signedSpeeds,frameRate,numSubSamples,sampleSizeFraction);
        plotbins = densitybins(2:end) - mean(diff(densitybins)); % convert from bin edges to centers
        revRates = mean(revRateSubsamples);
        revRateErr = std(revRateSubsamples);
        boundedline(plotbins,revRates,[revRateErr; revRateErr]',revFig.Children,'--')        
        xlabel(revFig.Children,['\rho_' num2str(k) ' (worms/mm^2)'])
        ylabel(revFig.Children,'relative reversal rate')
        figurename = ['figures/reversals/phaseSpecific/reversals-density_knn' num2str(k) ...
            '_pharynx_' strain '_' wormnum '_' phase '_data' num2str(dataset)];
        formatAndExportFigure(revFig,figurename,useJoinedTraj,exportOptions)
        %% plot histogram of increase in fwd frequency vs density
        [~, densitybins] = histcounts(knndensityAtFwd(:,k),'BinLimits',[0 densityBinLimits(strain)]);
        % subsample to estimate variability
        numSubSamples = 100;
        sampleSizeFraction = 0.05;
        fwdRateSubsamples = subsampleRevFreq(knndensityAtFwd(:,k),knndensity(:,k),...
            densitybins,-signedSpeeds,frameRate,numSubSamples,sampleSizeFraction);
        plotbins = densitybins(2:end) - mean(diff(densitybins)); % convert from bin edges to centers
        fwdRates = mean(fwdRateSubsamples);
        fwdRateErr = std(fwdRateSubsamples);
        boundedline(plotbins,fwdRates,[fwdRateErr; fwdRateErr]',fwdFig.Children,':')        
        xlabel(fwdFig.Children,['\rho_' num2str(k) ' (worms/mm^2)'])
        ylabel(fwdFig.Children,'relative forward rate')
        figurename = ['figures/reversals/phaseSpecific/forwards-density_knn' num2str(k) ...
            '_pharynx_' strain '_' wormnum '_' phase '_data' num2str(dataset)];
        formatAndExportFigure(fwdFig,figurename,useJoinedTraj,exportOptions)
        %% plot histogram of increase in rev or fwd frequency vs density
        knndensityEitherRevFwd = [knndensityAtRev(:,k); knndensityAtFwd(:,k)];
        [~, densitybins] = histcounts(knndensityEitherRevFwd,'BinLimits',[0 densityBinLimits(strain)]);
        % subsample to estimate variability
        numSubSamples = 100;
        sampleSizeFraction = 0.05;
        revfwdRateSubsamples = subsampleRevFreq(knndensityEitherRevFwd,knndensity(:,k),...
            densitybins,abs(signedSpeeds),frameRate,numSubSamples,sampleSizeFraction);
        plotbins = densitybins(2:end) - mean(diff(densitybins)); % convert from bin edges to centers
        revfwdRates = mean(revfwdRateSubsamples);
        revfwdRateErr = std(revfwdRateSubsamples);
        boundedline(plotbins,revfwdRates,[revfwdRateErr; revfwdRateErr]',revfwdFig.Children)        
        xlabel(revfwdFig.Children,['\rho_' num2str(k) ' (worms/mm^2)'])
        ylabel(revfwdFig.Children,'relative reorientation rate')
        figurename = ['figures/reversals/phaseSpecific/reorient-density_knn' num2str(k) ...
            '_pharynx_' strain '_' wormnum '_' phase '_data' num2str(dataset)];
        formatAndExportFigure(revfwdFig,figurename,useJoinedTraj,exportOptions)
        %% plot histogram of increase in reversal frequency vs change in density
        [~, densityChangebins] = histcounts(knndensityChangeAtRev(:,k),'BinLimits',[-2 2]*densityBinLimits(strain));
        % subsample to estimate variability
        numSubSamples = 100;
        sampleSizeFraction = 0.05;
        revRateSubsamples = subsampleRevFreq(knndensityChangeAtRev(:,k),knndensityChange(:,k),...
            densityChangebins,signedSpeeds,frameRate,numSubSamples,sampleSizeFraction);
        plotbins = densityChangebins(2:end) - mean(diff(densityChangebins)); % convert from bin edges to centers
        revRates = mean(revRateSubsamples);
        revRateErr = std(revRateSubsamples);
        boundedline(plotbins,revRates,[revRateErr; revRateErr]',revFig2.Children,'--')        
        xlabel(revFig2.Children,['\Delta_t\rho_' num2str(k) ' (worms/mm^2/s)'])
        ylabel(revFig2.Children,'relative reversal rate')
        figurename = ['figures/reversals/phaseSpecific/reversals-densityChange_knn' num2str(k) ...
            '_pharynx_' strain '_' wormnum '_' phase '_data' num2str(dataset)];
        formatAndExportFigure(revFig2,figurename,useJoinedTraj,exportOptions)
        %% plot histogram of increase in fwd frequency vs density
        [~, densityChangebins] = histcounts(knndensityChangeAtFwd(:,k),'BinLimits',[-2 2]*densityBinLimits(strain));
        % subsample to estimate variability
        numSubSamples = 100;
        sampleSizeFraction = 0.05;
        fwdRateSubsamples = subsampleRevFreq(knndensityChangeAtFwd(:,k),knndensityChange(:,k),...
            densityChangebins,-signedSpeeds,frameRate,numSubSamples,sampleSizeFraction);
        plotbins = densityChangebins(2:end) - mean(diff(densityChangebins)); % convert from bin edges to centers
        fwdRates = mean(fwdRateSubsamples);
        fwdRateErr = std(fwdRateSubsamples);
        boundedline(plotbins,fwdRates,[fwdRateErr; fwdRateErr]',fwdFig2.Children,':')        
        xlabel(fwdFig2.Children,['\Delta_t\rho_' num2str(k) ' (worms/mm^2/s)'])
        ylabel(fwdFig2.Children,'relative forward rate')
        figurename = ['figures/reversals/phaseSpecific/forwards-densityChange_knn' num2str(k) ...
            '_pharynx_' strain '_' wormnum '_' phase '_data' num2str(dataset)];
        formatAndExportFigure(fwdFig2,figurename,useJoinedTraj,exportOptions)
        %% plot histogram of increase in rev or fwd frequency vs density
        knndensityChangeEitherRevFwd = [knndensityChangeAtRev(:,k); knndensityChangeAtFwd(:,k)];
        [~, densityChangebins] = histcounts(knndensityChangeEitherRevFwd,'BinLimits',[-2 2]*densityBinLimits(strain));
        % subsample to estimate variability
        numSubSamples = 100;
        sampleSizeFraction = 0.05;
        revfwdRateSubsamples = subsampleRevFreq(knndensityChangeEitherRevFwd,knndensity(:,k),...
            densityChangebins,abs(signedSpeeds),frameRate,numSubSamples,sampleSizeFraction);
        plotbins = densityChangebins(2:end) - mean(diff(densityChangebins)); % convert from bin edges to centers
        revfwdRates = mean(revfwdRateSubsamples);
        revfwdRateErr = std(revfwdRateSubsamples);
        boundedline(plotbins,revfwdRates,[revfwdRateErr; revfwdRateErr]',revfwdFig2.Children)        
        xlabel(revfwdFig2.Children,['\Delta_t\rho_' num2str(k) ' (worms/mm^2/s)'])
        ylabel(revfwdFig2.Children,'relative reorientation rate')
        figurename = ['figures/reversals/phaseSpecific/reorient-densityChange_knn' num2str(k) ...
            '_pharynx_' strain '_' wormnum '_' phase '_data' num2str(dataset)];
        formatAndExportFigure(revfwdFig2,figurename,useJoinedTraj,exportOptions)
    end
    densityFig.Children.XLabel.String = ['\rho_' num2str(k) ' (worms/mm^2)'];
    densityFig.Children.XLabel.String = ['P'];
    lh = legend(densityFig.Children,num2str(knnbrNumValues'));
    lh.Title.String = 'k nearest nbrs';
    figurename = ['figures/reversals/phaseSpecific/density_convergence' ...
        '_pharynx_' strain '_' wormnum '_' phase '_data' num2str(dataset)];
    formatAndExportFigure(densityFig,figurename,useJoinedTraj,exportOptions)
end

end

function subsamples = subsampleRevFreq(knndensityAtRev,knndensity,densitybins,...
    signedSpeeds,frameRate,numSubSamples,sampleSizeFraction)
subsamples = NaN(numSubSamples,length(densitybins) - 1);
numSamplesRev = round(sampleSizeFraction*numel(knndensityAtRev));
numSamplesDensity = round(sampleSizeFraction*numel(knndensity));
for subSampleCtr = 1:numSubSamples
    sampleDensityAtRev = datasample(knndensityAtRev,numSamplesRev,'Replace',false);
    [sampleDensity, sampleIdcs] = datasample(knndensity,numSamplesDensity,'Replace',false);
    sampleNkAtRev = histcounts(sampleDensityAtRev,densitybins);
    sampleNkforward = histcounts(sampleDensity(signedSpeeds(sampleIdcs)>0),densitybins);
    sampleTkforward = sampleNkforward/frameRate; % convert frames to time
    sampleRevRate = sampleNkAtRev./sampleTkforward;
    sampleRevRate = sampleRevRate./sampleRevRate(1); % divide by lowest density estimate to get relative reversal rate
    subsamples(subSampleCtr,:) = sampleRevRate;
end
end

function [] = formatAndExportFigure(handle,figurename,useJoinedTraj,exportOptions)
set(handle,'PaperUnits','centimeters')
if useJoinedTraj, figurename = [figurename '_jointraj']; end
exportfig(handle,[figurename '.eps'],exportOptions)
system(['epstopdf ' figurename '.eps']);
system(['rm ' figurename '.eps']);
end