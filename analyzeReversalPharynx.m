function [] = analyzeReversalPharynx(dataset,phase,wormnum,useJoinedTraj)
% calculate various single worm statistics for different numbers of worms
% on plate
% INPUTS
% dataset: 1 or 2. To specify which dataset to run the script for.
% phase: 'joining', 'fullMovie', or 'sweeping'. Script defines stationary phase as: starts at 10% into the movie, and stops at 60% into the movie (HA and N2) or at specified stopping frames (npr-1).
% wormnum: '1W', '40', or 'HD'
% OUTPUTS
% none returned, but figures are exported

% issues / todo:

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);
load ~/Dropbox/Utilities/colormaps_ascii/increasing_cool/cmap_Blues.txt
addpath('auxiliary/')
addpath('filters/')
%% set parameters
if dataset ==1
    strains = {'N2','npr1'}; %{'npr1','HA','N2'}
elseif dataset ==2
    strains = {'N2','npr1'}; %{'npr1','N2'}
end
if dataset == 1
    intensityThresholds_g = containers.Map({'40','HD','1W'},{50, 40, 100});
elseif dataset ==2
    intensityThresholds_g = containers.Map({'40','HD','1W'},{60, 40, 100});
end
if nargin<4
    useJoinedTraj = true
end
maxBlobSize_g = 1e4;
minPathLength = 50; % minimum path length of reversals to be included - careful about choosing this too high, as it will also apply to incompletely tracked reversals
knnbrNumValues = 6; % which numbers of knnbrs to consider for checking density dependence
numk = numel(knnbrNumValues);
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
minSpeedPerFrame = pixelsize/2; % minimum speed to consider for reversals, per frame. worms slower than this are considered paused and not moving forward or backward
speedBinLimits = containers.Map({'npr1','N2'},{400, 300});
densityBinLimits = containers.Map({'npr1','N2'},{10.5, 5.5});
speedThreshold = containers.Map({'npr1','N2'},{100,50});
nStrains = length(strains);
plotColors = flipud(lines(nStrains));
% lineHandles = NaN(nStrains,1);
densityBinWidth = 0.25;
nDensityBins = 10; % for plotting speed histograms
speedFractions = NaN(nStrains,2*nDensityBins,2);
speedBinWidth = 10;
numSubSamples = 100;
sampleSizeFraction = 1;
samplingCorrectionFactor = sqrt(sampleSizeFraction); % sqrt(b/n), see Geyer 2013, subsampling bootstrap
if sampleSizeFraction < 1
    sampleReplacement = false;
else
    sampleReplacement = true;
end

%% go through strains, densities, movies
for strainCtr = 1:nStrains
    strain = strains{strainCtr};
    %% load data
    if dataset ==1
        [phaseFrames,filenames,~] = xlsread(['datalists/' strain '_' wormnum '_list_lslx.xlsx'],1,'A1:E15','basic');
    elseif dataset ==2
        [phaseFrames,filenames,~] = xlsread(['datalists/' strain '_' wormnum '_g_list_lslx.xlsx'],1,'A1:E15','basic');
    end
    %     if ~useJoinedTraj
    %         filenames = strrep(filenames,'/data2/shared/data/twoColour/Results/',...
    %             '/data2/databackup/data/twoColour/ResultsUnjoinedTrajectories/');
    %     end
    %% intialize variables
    numFiles = length(filenames);
    knndensity = cell(numFiles,1);
    knndensityAtRev = cell(numFiles,1);
    knndensityAtFwd = cell(numFiles,1);
    signedSpeeds = cell(numFiles,1);
    %% loop through files and load data
    %%%%
    numFramesAllFiles = 0;
    %%%%%
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
            %%%
            numFramesAllFiles = numFramesAllFiles + nnz(nonNanFramesLogInd);
            disp(['total frame count ' num2str(numFramesAllFiles)])
            %%%%
            signedSpeeds{fileCtr} = signedSpeedThisFile(nonNanFramesLogInd);
            %% find reversals in signed speed
            [revStartInd, ~, untrackedRevEnds, ~, incompleteInterRev] = ...
                findReversals(signedSpeedThisFile,trajData.worm_index_joined,minPathLength,frameRate,minSpeedPerFrame);
            [fwdStartInd, ~, ~, ~, ~] = ...
                findReversals(-signedSpeedThisFile,trajData.worm_index_joined,minPathLength,frameRate,minSpeedPerFrame);
            disp([num2str(100*mean(incompleteInterRev),2) '% of tracks are lost between reversals'])
            disp([num2str(100*mean(untrackedRevEnds),2) '% of tracks are lost during reversals'])
            % also set all reversals with untracked ends as incomplete interrevs
            incompleteInterRev = incompleteInterRev|untrackedRevEnds;
            disp([num2str(100*mean(incompleteInterRev),2) '% of reversals are thus not fully tracked' newline])
            %% estimate density based on k-nearest neighbours
            neighbr_dist = h5read(filename,'/neighbr_distances')/1000;% convert to mm
            knndensityThisFile = NaN(numFramesTotal,10);
            for k = 1:10
                knndensityThisFile(:,k) = k./(pi*neighbr_dist(:,k).^2);
                % frames with zero density result from the nbr-distance being
                % Inf, less than k worms tracked for that frame. Hence, we want
                % to discard these (an alternative would be to use the non-zero
                % estimate from a lower k)
                knndensityThisFile(isinf(neighbr_dist(:,k)),k) = NaN;
            end
            % ignore frames where worm index changes
            knndensityThisFile(wormChangeIndcs,:) = NaN;
            % smooth density for 1s
            knndensityThisFile = smoothdata(knndensityThisFile,1,'movmean',round(frameRate),'omitnan');
            knndensity{fileCtr} = knndensityThisFile(nonNanFramesLogInd,:);
            %% store reversal statistics for this file
            % densities at reversal
            knndensityAtRev{fileCtr} = knndensityThisFile(revStartInd,:);
            knndensityAtFwd{fileCtr} = knndensityThisFile(fwdStartInd,:);
        end
    end
    %% pool data from all files
    knndensity = vertcat(knndensity{:});
    knndensityAtRev = vertcat(knndensityAtRev{:});
    knndensityAtFwd = vertcat(knndensityAtFwd{:});
    signedSpeeds = vertcat(signedSpeeds{:});
    %% loop over density estimates
    %     densityFig = figure; hold on
    speedAbsHist2Dfig = figure; hold on
    colormap(speedAbsHist2Dfig,flipud(cmap_Blues))
    for k = knnbrNumValues
        % initialize figures
        if strainCtr==1
            revFig = figure; hold on
            fwdFig = figure; hold on
            revfwdFig = figure; hold on
        end
        %         %% plot "convergence" of knn density estimate
        %         histogram(densityFig.Children, knndensity(:,k),'Normalization','Probability','DisplayStyle','stairs','BinLimits',[0 20])
        %% plot 2d histogram of absolute speed vs density
        h = histogram2(speedAbsHist2Dfig.Children,knndensity(:,k),abs(signedSpeeds),...
            'DisplayStyle','tile','YBinLimits',[minSpeedPerFrame*frameRate speedBinLimits(strain)],...
            'BinWidth',[densityBinWidth, speedBinWidth],'XBinLimits',[0 densityBinLimits(strain)],'EdgeColor','none','Normalization','Probability');
        normfactor = sum(h.BinCounts,2); % for conditional normalisation
        normfactor(normfactor==0) = 1;
        h.BinCounts = h.BinCounts./normfactor; % conditional normalisation
        xlabel(speedAbsHist2Dfig.Children,['local density, \rho_' num2str(k) ' (worms/mm^2)'])
        ylabel(speedAbsHist2Dfig.Children,'speed (\mum/s)')
        figurename = ['figures/reversals/phaseSpecific/speedabs-density_knn' num2str(k) ...
            '_pharynx_' strain '_' wormnum '_' phase '_data' num2str(dataset)];
        %         formatAndExportFigure(speedAbsHist2Dfig,figurename,useJoinedTraj,exportOptions)
        %% plot speed histograms for first few density bins
        if k==6
            speedAbsHistFig = figure; hold on
            speedAbsHistFig.Children.ColorOrder = flipud(cool(nDensityBins));
            for binCtr = 1:(2*nDensityBins)
                thisBinLogInd = knndensity(:,k)>=h.XBinEdges(binCtr)&knndensity(:,k)<h.XBinEdges(binCtr+1);
                [thisbincounts, thisbinedges] = histcounts(abs(signedSpeeds(thisBinLogInd)),...
                    'BinWidth',2*speedBinWidth,'Normalization','Probability');
                if binCtr<=nDensityBins
                    plot(speedAbsHistFig.Children,thisbinedges(1:end-1)+speedBinWidth,thisbincounts);
                end
                % calculate proportion in high and low speed states
                speedFractions(strainCtr,binCtr,1) = sum(thisbincounts(thisbinedges(2:end)<=speedThreshold(strain)));
                speedFractions(strainCtr,binCtr,2) = sum(thisbincounts(thisbinedges(1:end-1)>speedThreshold(strain)...
                    &thisbinedges(2:end)<=350));
            end
            %             speedAbsHistFig.Children.Box = 'on';
            %             if strcmp(strain,'npr1')
            %                 speedAbsHistFig.Children.XLim(2) = 600;
            %             elseif strcmp(strain,'N2')
            %                 speedAbsHistFig.Children.XLim(2) = 400;
            %             end
            %             speedAbsHistFig.Children.YTick = [0 0.05 0.1 0.15 0.2];
            %             xlabel(speedAbsHistFig.Children,'speed (\mum/s)')
            %             ylabel(speedAbsHistFig.Children,'relative frequency')
            %             % add a colorbar
            %             hcb = colorbar;
            %             hcb.Location = 'east';
            %             colormap(flipud(cool(nDensityBins)));
            %             hcb.Ticks = [0 1];
            %             hcb.Label.String = 'local density';
            %             hcb.TickLabels = {'0',num2str(h.XBinEdges(nDensityBins+1))};
            %             hcb.Position(4) = 0.5;
            %             hcb.Position(2) = 0.25;
            %             % export figure
            %             figurename = ['figures/reversals/phaseSpecific/speedabs-density_1D_knn' num2str(k) ...
            %                 '_pharynx_' strain '_' wormnum '_' phase '_data' num2str(dataset)];
            %             formatAndExportFigure(speedAbsHistFig,figurename,useJoinedTraj,exportOptions)
        end
        %% plot increase in reversal frequency vs density
        [~, densitybins] = histcounts(knndensityAtRev(:,k),'BinWidth',2*densityBinWidth,'BinLimits',[0 densityBinLimits(strain)]);
        % subsample to estimate variability
        revRateSubsamples = subsampleRevFreq(knndensityAtRev(:,k),knndensity(:,k),...
            densitybins,signedSpeeds,frameRate,minSpeedPerFrame,numSubSamples,sampleSizeFraction,sampleReplacement);
        %%% careful: last line uses frameRate of last files, only works if all frameRates are consistent
        plotbins = densitybins(2:end) - diff(densitybins)./2; % convert from bin edges to centers
        revRates = mean(revRateSubsamples,1);
        revRateErr = std(revRateSubsamples,0,1).*samplingCorrectionFactor;
        [~, mindensityIdx] = min(abs(densitybins));
        revRatesNormalised = revRates./(revRates(mindensityIdx) + eps); % divide by lowest density estimate to get relative reversal rate
        revRateErrNormalised = revRatesNormalised.*sqrt((revRateErr./(revRates + eps)).^2 + ...
            (revRateErr(mindensityIdx)./(revRates(mindensityIdx) + eps)).^2); % sum of square relative errors rule
        % this normalisation could result in large estimates for low outliers
        % (where revRates(mindensityIdx)==0), so may it be better to
        % normalise by the mean, or the lowest non-zero estimate?
        boundedline(plotbins,revRatesNormalised,[revRateErrNormalised; revRateErrNormalised]',...
            revFig.Children,'-','cmap',plotColors(strainCtr,:),'alpha')
        %                 boundedline(plotbins,revRates,[revRateErr; revRateErr]',...
        %             revFig.Children,'-','cmap',plotColors(strainCtr,:),'alpha')
        if strainCtr==nStrains
            xlabel(revFig.Children,['local density, \rho_' num2str(k) ' (worms/mm^2)'])
            ylabel(revFig.Children,'relative reversal rate')
            revFig.Children.YLim(1) = 0;
            revFig.Children.XLim(2) = 5;
            revFig.Children.Box = 'on';
%             % adjust aspect ratio
%             revFig.Position(4) = 3/4*revFig.Position(3);
            figurename = ['figures/reversals/phaseSpecific/reversals-density_knn' num2str(k) ...
                '_pharynx_' wormnum '_' phase '_data' num2str(dataset) '_minPath' num2str(minPathLength)];
            formatAndExportFigure(revFig,figurename,useJoinedTraj,exportOptions)
            % now make second version with inferred revrate slopes for reference
            % npr-1 inferred revrate
            plot(revFig.Children,plotbins,revRatesNormalised(mindensityIdx) ...
                + (plotbins-densityBinWidth/2)*0.265*0.6234,'--','Color',plotColors(2,:),'LineWidth',2)
            % N2 inferred revrate
            plot(revFig.Children,plotbins,revRatesNormalised(mindensityIdx)...
                + (plotbins-densityBinWidth/2)*0.064*0.6234,'--','Color',plotColors(1,:),'LineWidth',2)
            figurename = ['figures/reversals/phaseSpecific/reversals-density_knn' num2str(k) ...
                '_pharynx_' wormnum '_' phase '_data' num2str(dataset) '_minPath' num2str(minPathLength) '_modelComparison'];
            %         legend({strain,'$r_\mathrm{rev}(\rho)$'},'Location','NorthWest','Interpreter','Latex')
%             % adjust aspect ratio
%             revFig.Position(4) = 3/4*revFig.Position(3);
            formatAndExportFigure(revFig,figurename,useJoinedTraj,exportOptions)
        end
        
        %         %% plot increase in fwd frequency vs density
        %         [~, densitybins] = histcounts(knndensityAtFwd(:,k),'BinWidth',densityBinWidth,'BinLimits',[0 densityBinLimits(strain)]);
        %         % subsample to estimate variability
        %         fwdRateSubsamples = subsampleRevFreq(knndensityAtFwd(:,k),knndensity(:,k),...
        %             densitybins,-signedSpeeds,frameRate,minSpeedPerFrame,numSubSamples,sampleSizeFraction,sampleReplacement);
        %         %%% careful: last line uses frameRate of last files, only works if all frameRates are consistent
        %         plotbins = densitybins(2:end) - diff(densitybins)./2; % convert from bin edges to centers
        %         fwdRates = mean(fwdRateSubsamples,1);
        %         fwdRateErr = std(fwdRateSubsamples,0,1).*samplingCorrectionFactor;
        %         [~, mindensityIdx] = min(abs(densitybins));
        %         fwdRatesNormalised = fwdRates./(fwdRates(mindensityIdx) + eps); % divide by lowest density estimate to get relative reversal rate
        %         fwdRateErrNormalised = fwdRatesNormalised.*sqrt((fwdRateErr./(fwdRates + eps)).^2 + ...
        %             (fwdRateErr(mindensityIdx)./(fwdRates(mindensityIdx) + eps)).^2); % sum of square relative errors rule
        %         boundedline(plotbins,fwdRatesNormalised,[fwdRateErrNormalised; fwdRateErrNormalised]',...
        %             fwdFig.Children,':','cmap',plotColors(strainCtr,:),'alpha')
        %         if strainCtr==nStrains
        %             xlabel(fwdFig.Children,['local density, \rho_' num2str(k) ' (worms/mm^2)'])
        %             ylabel(fwdFig.Children,'relative forward rate')
        %             fwdFig.Children.YLim(1) = 0;
        %             fwdFig.Children.XLim(2) = 5;
        %             fwdFig.Children.Box = 'on';
        %             figurename = ['figures/reversals/phaseSpecific/forwards-density_knn' num2str(k) ...
        %                 '_pharynx_' wormnum '_' phase '_data' num2str(dataset) '_minPath' num2str(minPathLength)];
        %             formatAndExportFigure(fwdFig,figurename,useJoinedTraj,exportOptions)
        %         end
        %         %% plot increase in rev or fwd frequency vs density
        %         knndensityEitherRevFwd = [knndensityAtRev(:,k); knndensityAtFwd(:,k)];
        %         [~, densitybins] = histcounts(knndensityEitherRevFwd,'BinWidth',densityBinWidth,'BinLimits',[0 densityBinLimits(strain)]);
        %         % subsample to estimate variability
        %         revfwdRateSubsamples = subsampleRevFreq(knndensityEitherRevFwd,knndensity(:,k),...
        %             densitybins,abs(signedSpeeds),frameRate,minSpeedPerFrame,numSubSamples,sampleSizeFraction,sampleReplacement);
        %         %%% careful: last line uses frameRate of last files, only works if all frameRates are consistent
        %         plotbins = densitybins(2:end) - diff(densitybins)./2; % convert from bin edges to centers
        %         revfwdRates = mean(revfwdRateSubsamples,1);
        %         revfwdRateErr = std(revfwdRateSubsamples,0,1).*samplingCorrectionFactor;
        %         [~, mindensityIdx] = min(abs(densitybins));
        %         revfwdRatesNormalised = revfwdRates./(revfwdRates(mindensityIdx) + eps); % divide by lowest density estimate to get relative reversal rate
        %         revfwdRateErrNormalised = revfwdRatesNormalised.*sqrt((revfwdRateErr./(revfwdRates + eps)).^2 + ...
        %             (revfwdRateErr(mindensityIdx)./(revfwdRates(mindensityIdx) + eps)).^2); % sum of square relative errors rule
        %         boundedline(plotbins,revfwdRatesNormalised,[revfwdRateErrNormalised; revfwdRateErrNormalised]',...
        %             revfwdFig.Children,'cmap',plotColors(strainCtr,:),'alpha')
        %         if strainCtr==nStrains
        %             xlabel(revfwdFig.Children,['local density, \rho_' num2str(k) ' (worms/mm^2)'])
        %             ylabel(revfwdFig.Children,'relative reorientation rate')
        %             revfwdFig.Children.YLim(1) = 0;
        %             revfwdFig.Children.XLim(2) = 5;
        %             revfwdFig.Children.Box = 'on';
        %             figurename = ['figures/reversals/phaseSpecific/reorient-density_knn' num2str(k) ...
        %                 '_pharynx_' wormnum '_' phase '_data' num2str(dataset) '_minPath' num2str(minPathLength)];
        %             formatAndExportFigure(revfwdFig,figurename,useJoinedTraj,exportOptions)
        %         end
        %
    end
    %     densityFig.Children.XLabel.String = ['local density, \rho_' num2str(k) ' (worms/mm^2)'];
    %     densityFig.Children.XLabel.String = ['P'];
    %     lh = legend(densityFig.Children,num2str(knnbrNumValues'));
    %     lh.Title.String = 'k nearest nbrs';
    %     figurename = ['figures/reversals/phaseSpecific/density_convergence' ...
    %         '_pharynx_' strain '_' wormnum '_' phase '_data' num2str(dataset)];
    %     formatAndExportFigure(densityFig,figurename,useJoinedTraj,exportOptions)
end

% plot ratio of speed fractions vs density
figure
thisDensityBins = 0:densityBinWidth:(2*nDensityBins-1)*densityBinWidth;
plot(thisDensityBins,speedFractions(1,:,2)./speedFractions(1,:,1),'Color',plotColors(1,:))
hold on
plot(thisDensityBins,speedFractions(2,:,2)./speedFractions(2,:,1),'Color',plotColors(2,:))
ylabel('P_f/P_s')
xlabel(['local density, \rho_' num2str(k) ' (worms/mm^2)'])
legend(strains,'Location','NorthEast')
formatAndExportFigure(gcf,'speedRatios',false,exportOptions)
end

function subsamples = subsampleRevFreq(knndensityAtRev,knndensity,densitybins,...
    signedSpeeds,frameRate,minSpeedPerFrame,numSubSamples,sampleSizeFraction,replacement)
subsamples = NaN(numSubSamples,length(densitybins) - 1);
numSamplesRev = round(sampleSizeFraction*numel(knndensityAtRev));
numSamplesDensity = round(sampleSizeFraction*numel(knndensity));
for subSampleCtr = 1:numSubSamples
    % count reversal events for a given density
    sampleDensityAtRev = datasample(knndensityAtRev,numSamplesRev,'Replace',replacement);
    % compare this to number of frames at the same densities
    [sampleDensity, sampleIdcs] = datasample(knndensity,numSamplesDensity,'Replace',replacement);
    sampleNkAtRev = histcounts(sampleDensityAtRev,densitybins);
    % restrict to those frames with the right movement direction (eg forward)
    sampleNkforward = histcounts(sampleDensity(signedSpeeds(sampleIdcs)>minSpeedPerFrame*frameRate),densitybins);
    sampleTkforward = sampleNkforward/frameRate; % convert frames to time
    % estimate reversal rate as number of events per time (Poisson rate estimate)
    sampleRevRate = sampleNkAtRev./(sampleTkforward);
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