function [] = analyzeCorrelationData(dataset,phase,wormnum,markerType,plotDiagnostics)
% calculate speed vs neighbr distance, directional correlation, and
% radial distribution functions
% INPUTS
% dataset: 1 or 2. To specify which dataset to run the script for.
% phase: 'joining', 'fullMovie', or 'sweeping'. Script defines stationary phase as: starts at 10% into the movie, and stops at 60% into the movie (HA and N2) or at specified stopping frames (npr-1).
% wormnum: '40', or 'HD'
% markerType: 'pharynx', or 'bodywall'
% plotDiagnostics: true (default) or false
% OUTPUTS
% none returned, but figures are exported
% issues/to-do:
% - seperate into individual functions for each statistic?
% - calculate red-green correlations as well as red-red

addpath('auxiliary/')
addpath('filters/')

%% set fixed parameters

if nargin<5
    plotDiagnostics = false; % true or false
    if nargin<4
        markerType = 'pharynx';
    end
end

if dataset ==1
    strains = {'npr1','N2'};%{'npr1','HA','N2'}
    assert(~strcmp(markerType,'bodywall'),'Bodywall marker for dataset 1 not available')
elseif dataset ==2
    strains = {'npr1','N2'};
end
useJoinedTraj = false;

nStrains = length(strains);

% filtering parameters
if dataset == 1
    intensityThresholds = containers.Map({'40','HD','1W'},{50, 40, 100});
elseif dataset ==2
    intensityThresholds = containers.Map({'40','HD','1W'},{60, 40, 100});
end
if strcmp(markerType,'pharynx')
    maxBlobSize = 1e5;
    channelStr = 'g';
elseif strcmp(markerType,'bodywall')
    maxBlobSize = 2.5e5;
    channelStr = 'r';
    minSkelLength = 850;
    maxSkelLength = 1500;
else
    error('unknown marker type specified, should be pharynx or bodywall')
end
% analysis parameters
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
distBinWidth = 100; % in units of micrometers
maxSpeed = 500;
minSpeedPerFrame = pixelsize/2; % assume these worms to non-moving, and exclude from velocity correlation analysis (as that will be inaccurate for very low magnitude velocities)
maxDist = 2000;
minDist = 70;
distBins = 0:distBinWidth:maxDist;

% plotting parameters
plotColors = lines(nStrains);
if plotDiagnostics, visitfreqFig = figure; hold on, end
distxticks = 0:0.5:(maxDist/1000);
load('~/Dropbox/Utilities/colormaps_ascii/increasing_cool/cmap_Blues.txt')
% export fig parameters
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);
%% initialize figures
dircorrFig = figure; hold on
velcorrFig = figure; hold on
velnbrcorrFig = figure; hold on
absvelnbrcorrFig = figure; hold on
accnbrcorrFig = figure; hold on
poscorrFig = figure; hold on
orderFig = figure; hold on
lineHandles = NaN(nStrains,1);
%% loop through strains
for strainCtr = 1:nStrains
    angleFig = figure; hold on
    %% load file lists
    if dataset == 1
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_list_lslx.xlsx'],1,'A1:E15','basic');
    elseif dataset == 2
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_' channelStr '_list_lslx.xlsx'],1,'A1:E15','basic');
    end
%     if ~useJoinedTraj
%         filenames = strrep(filenames,'/data2/shared/data/twoColour/Results/',...
%             '/end/home/lschumac/databackup/data/twoColour/ResultsUnjoinedTrajectories/');
%     end
    numFiles = length(filenames);
    if strcmp(wormnum,'40'), visitfreq = cell(numFiles,1); end
    %% intialise variables to store results
    dirxcorr = cell(numFiles,1); % for calculating directional cross-correlation
    velxcorr = cell(numFiles,1); % for calculating velocity cross-correlation
    velnbrcorr = cell(numFiles,1); % for calculating velocity correlation with nearest neighbor position
    velnbrcorrSigned = cell(numFiles,1); % for calculating velocity correlation with nearest neighbor position
    velnbrcorrFwd = cell(numFiles,1); % for calculating velocity correlation with nearest neighbor position
    velnbrcorrCoM = cell(numFiles,1); % for calculating velocity correlation with nearest neighbor position
    accnbrcorr = cell(numFiles,1);
    pairdist = cell(numFiles,1);
    nbrDist= cell(numFiles,1); 
    distCoM = cell(numFiles,1);
    polar_order= cell(numFiles,1);
    nematic_order= cell(numFiles,1);
    pcf =cell(numFiles,1);
    %% loop through files
    for fileCtr = 1:numFiles % can be parfor
        filename = filenames{fileCtr};
        %% load tracking data
        trajData = h5read(filename,'/trajectories_data');
        blobFeats = h5read(filename,'/blob_features');
        skelData = h5read(filename,'/skeleton');
        % check formats
        assert(size(skelData,1)==2,['Wrong skeleton size for ' filename])
        assert(size(skelData,3)==length(trajData.frame_number),['Wrong number of skeleton frames for ' filename])
        if all(isnan(skelData(:))), warning(['all skeleton are NaN for ' filename]),end
        %% randomly sample which frames to analyze
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr);
        numFrames = round((lastFrame-firstFrame)/frameRate/3);
        framesAnalyzed = randperm((lastFrame-firstFrame),numFrames) + firstFrame; % randomly sample frames without replacement
        %% filter data for worms
        minSpeed = minSpeedPerFrame*frameRate;
        if plotDiagnostics
            visualizeIntensitySizeFilter(blobFeats,pixelsize,intensityThresholds(wormnum),maxBlobSize,...
                [wormnum ' ' strains{strainCtr} ' ' strrep(filename(end-32:end-18),'/','')])
        end
        if strcmp(markerType,'pharynx')
            % reset skeleton flag for pharynx data
            trajData.has_skeleton = true(size(trajData.has_skeleton)); %squeeze(~any(any(isnan(skelData))));
        end
        trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
            intensityThresholds(wormnum),maxBlobSize)...
            &trajData.has_skeleton; % careful: we may not want to filter for skeletonization for
        %clustering statistics
        if strcmp(markerType,'bodywall')
            % filter red data by skeleton length
            trajData.filtered = trajData.filtered&logical(trajData.is_good_skel)&...
                filterSkelLength(skelData,pixelsize,minSkelLength,maxSkelLength);
        end
        % apply phase restriction
        phaseFilter_logInd = trajData.frame_number < lastFrame & trajData.frame_number > firstFrame;
        trajData.filtered(~phaseFilter_logInd)=false;
        %% calculate stats
        % calculate area for normalisation of pair correlation function
        if strcmp(wormnum,'40')
            OverallArea = pi*(8500/2)^2;
        else
            OverallArea = peak2peak(trajData.coord_x(trajData.filtered)).*...
                peak2peak(trajData.coord_y(trajData.filtered)).*pixelsize.^2;
            disp(['overall background area estimated as ' num2str(OverallArea)])
        end
        % initialise variables to store results for this file
        dirxcorr{fileCtr} = cell(numFrames,1); % for calculating directional cross-correlation
        velxcorr{fileCtr} = cell(numFrames,1); % for calculating velocity cross-correlation
        velnbrcorr{fileCtr} = cell(numFrames,1);
        velnbrcorrSigned{fileCtr} = cell(numFrames,1);
        velnbrcorrFwd{fileCtr} = cell(numFrames,1);
        velnbrcorrCoM{fileCtr} = cell(numFrames,1);
        accnbrcorr{fileCtr} = cell(numFrames,1);
        pairdist{fileCtr} = cell(numFrames,1);
        nbrDist{fileCtr}= cell(numFrames,1);
        distCoM{fileCtr}= cell(numFrames,1);
        pcf{fileCtr} = NaN(length(distBins) - 1,numFrames);
        polar_order{fileCtr} = zeros(numFrames,1);
        nematic_order{fileCtr} = zeros(numFrames,1);
        if strcmp(markerType,'bodywall')
            [ ~, velocities_x, velocities_y, ~, acceleration_x, acceleration_y ] = ...
                calculateSpeedsFromSkeleton(trajData,skelData,1:5,...
                pixelsize,frameRate,true,0);
        end
        for frameCtr = 1:numFrames % one may be able to vectorise this
            frame = framesAnalyzed(frameCtr);
            [x ,y] = getWormPositions(trajData, frame, true);
            x = double(x); y = double(y);
            N = length(x);
            if N>1 % need at least two worms in frame
                frameLogInd = trajData.frame_number==frame&trajData.filtered;
                if strcmp(markerType,'pharynx')
                    vx = double(blobFeats.velocity_x(frameLogInd))*pixelsize*frameRate;
                    vy = double(blobFeats.velocity_y(frameLogInd))*pixelsize*frameRate;
                    accx = double(blobFeats.acceleration_x(frameLogInd)); % units may need correcting
                    accy = double(blobFeats.acceleration_y(frameLogInd));
                    orientation_x = double(squeeze(skelData(1,1,frameLogInd) - skelData(1,2,frameLogInd)));
                    orientation_y = double(squeeze(skelData(2,1,frameLogInd) - skelData(2,2,frameLogInd)));
                elseif strcmp(markerType,'bodywall')
                    vx = velocities_x(frameLogInd)*pixelsize*frameRate;
                    vy = velocities_y(frameLogInd)*pixelsize*frameRate;
                    accx = acceleration_x(frameLogInd); % units may need correcting
                    accy = acceleration_y(frameLogInd);
                    orientation_x = double(squeeze(skelData(1,1,frameLogInd) - skelData(1,8,frameLogInd)));
                    orientation_y = double(squeeze(skelData(2,1,frameLogInd) - skelData(2,8,frameLogInd)));
                end
                % normalise orientation for head/pharynx size
                headSize = sqrt(orientation_x.^2 + orientation_y.^2);
                orientation_x = orientation_x./headSize;
                orientation_y = orientation_y./headSize;
                % assume speeds>maxSpeed are tracking errors, and
                % speeds<minSpeed to be below reasonable detectabilty
                % threshold
                v = sqrt(vx.^2 + vy.^2);
                speedOutlierLogIdcs = v>maxSpeed|v<minSpeed;
                if any(speedOutlierLogIdcs)
                    vx(speedOutlierLogIdcs) = NaN;
                    vy(speedOutlierLogIdcs) = NaN;
                    accx(speedOutlierLogIdcs) = NaN;
                    accy(speedOutlierLogIdcs) = NaN;
                end
                % check when worm is moving forward vs backward
                moveStateThisFrame = sign(smoothdata(double(blobFeats.signed_speed(frameLogInd)),...
                    'movmean',round(frameRate/2),'omitnan')); % smooth speed to denoise

                %% calculate dir-dir and vel-vel correlations
                dirxcorr{fileCtr}{frameCtr} = vectorCrossCorrelation2D(orientation_x,orientation_y,false,false); % directional correlation
                velxcorr{fileCtr}{frameCtr} = vectorCrossCorrelation2D(vx,vy,true,false); % velocity correlation
                %% calculate pairwise distances
                pairdist{fileCtr}{frameCtr} = pdist([x y]).*pixelsize; % distance between all pairs, in micrometer
                %% calculate pair correlation
                pcf{fileCtr}(:,frameCtr) = histcounts(pairdist{fileCtr}{frameCtr},distBins,'Normalization','count'); % radial distribution function
                pcf{fileCtr}(:,frameCtr) = pcf{fileCtr}(:,frameCtr)'.*OverallArea ...
                    ./(pi*(distBins(2:end).^2 - (distBins(2:end) - distBinWidth).^2)*N*(N-1)/2); % normalisation by N(N-1)/2 as pdist doesn't double-count pairs
                % keep only pair-distances below maxDist
                pairDistKeepLogInd = pairdist{fileCtr}{frameCtr}<=maxDist;
                pairdist{fileCtr}{frameCtr} = pairdist{fileCtr}{frameCtr}(pairDistKeepLogInd);
                dirxcorr{fileCtr}{frameCtr} = dirxcorr{fileCtr}{frameCtr}(pairDistKeepLogInd);
                velxcorr{fileCtr}{frameCtr} = velxcorr{fileCtr}{frameCtr}(pairDistKeepLogInd);
                %% calculate direction towards neighbor(s) for each worm
                dx = (x - x')*pixelsize; dy = (y - y')*pixelsize;
                dD = sqrt(dx.^2 + dy.^2);
                % keep only nbr-distances below maxDist
                nbrKeepLogIdcs = dD(:)<=maxDist;
                
                nbrDist{fileCtr}{frameCtr} = dD(nbrKeepLogIdcs);
                velnbrcorr{fileCtr}{frameCtr} = vectorPairedCorrelation2D(vx,vy,dx./dD,dy./dD,true,false);

                velnbrcorrSigned{fileCtr}{frameCtr} = velnbrcorr{fileCtr}{frameCtr}.*moveStateThisFrame;
                velnbrcorrSigned{fileCtr}{frameCtr}(moveStateThisFrame==0,:) = NaN;

                velnbrcorrFwd{fileCtr}{frameCtr} = velnbrcorr{fileCtr}{frameCtr};
                velnbrcorrFwd{fileCtr}{frameCtr}(moveStateThisFrame<=0,:) = NaN;
                
                velnbrcorrCoM{fileCtr}{frameCtr} = vectorPairedCorrelation2D(vx,vy,x - nanmean(x),y - nanmean(y),true,false);
                distCoM{fileCtr}{frameCtr} = sqrt((x - nanmean(x)).^2 + (y - nanmean(y)).^2);
                
                % keep only nbrs < maxDist and reshape to column vector
                velnbrcorr{fileCtr}{frameCtr} = velnbrcorr{fileCtr}{frameCtr}(nbrKeepLogIdcs);
                velnbrcorrSigned{fileCtr}{frameCtr} = velnbrcorrSigned{fileCtr}{frameCtr}(nbrKeepLogIdcs);
                velnbrcorrFwd{fileCtr}{frameCtr} = velnbrcorrFwd{fileCtr}{frameCtr}(nbrKeepLogIdcs);

                % calculate correlation between change in direction and nbr positions
                accnbrcorr{fileCtr}{frameCtr} = double(vectorPairedCorrelation2D(accx,accy,dx./dD,dy./dD,true,false));
                accnbrcorr{fileCtr}{frameCtr} = accnbrcorr{fileCtr}{frameCtr}(nbrKeepLogIdcs);
                %% check consistency of outputs
                if (numel(velnbrcorr{fileCtr}{frameCtr})~=numel(nbrDist{fileCtr}{frameCtr}))||...
                        (numel(dirxcorr{fileCtr}{frameCtr})~=numel(pairdist{fileCtr}{frameCtr}))
                    error(['Inconsistent number of variables in frame ' num2str(frame) ' of ' filename ])
                end
                %% calculate polar and nematic global order
                phis = atan2(orientation_y,orientation_x);
                polar_order{fileCtr}(frameCtr) = abs(nanmean(exp(1i*phis)));
                nematic_order{fileCtr}(frameCtr) = abs(nanmean(exp(2i*phis)));
            end
        end
        %% pool data from frames
        dirxcorr{fileCtr} = horzcat(dirxcorr{fileCtr}{:});
        velxcorr{fileCtr} = horzcat(velxcorr{fileCtr}{:});
        velnbrcorr{fileCtr} = vertcat(velnbrcorr{fileCtr}{:});
        velnbrcorrSigned{fileCtr} = vertcat(velnbrcorrSigned{fileCtr}{:});
        velnbrcorrFwd{fileCtr} = vertcat(velnbrcorrFwd{fileCtr}{:});
        velnbrcorrCoM{fileCtr} = vertcat(velnbrcorrCoM{fileCtr}{:});
        accnbrcorr{fileCtr} = vertcat(accnbrcorr{fileCtr}{:});
        pairdist{fileCtr} = horzcat(pairdist{fileCtr}{:});
        nbrDist{fileCtr} = vertcat(nbrDist{fileCtr}{:});
        distCoM{fileCtr} = vertcat(distCoM{fileCtr}{:});
        %% heat map of sites visited - this only makes sense for 40 worm dataset where we don't move the camera
        if strcmp(wormnum,'40')&& plotDiagnostics
            siteVisitFig = figure;
            h=histogram2(trajData.coord_x*pixelsize/1000,trajData.coord_y*pixelsize/1000,...
                'DisplayStyle','tile','EdgeColor','none','Normalization','pdf');
            visitfreq{fileCtr} = h.Values(:);
            cb = colorbar; cb.Label.String = '# visited';
            axis equal
            xlabel('x (mm)'), ylabel('y (mm)')
            title([strains{strainCtr} ' ' strrep(filename(end-32:end-18),'/','')])
            set(siteVisitFig,'PaperUnits','centimeters')
            figurename = ['figures/individualRecordings/' strains{strainCtr}...
                '_' strrep(strrep(filename(end-32:end-18),' ',''),'/','') '_sitesVisited' '_' phase '_data' num2str(dataset)];
            exportfig(siteVisitFig,[figurename '.eps'],exportOptions)
            system(['epstopdf ' figurename '.eps']);
            system(['rm ' figurename '.eps']);
        end
    end
    %% combine data from multiple files
    nbrDist = vertcat(nbrDist{:});
    distCoM = vertcat(distCoM{:});
    pairdist = horzcat(pairdist{:});
    dirxcorr = horzcat(dirxcorr{:});
    velxcorr = horzcat(velxcorr{:});
    velnbrcorr = vertcat(velnbrcorr{:});
    velnbrcorrSigned = vertcat(velnbrcorrSigned{:});
    velnbrcorrFwd = vertcat(velnbrcorrFwd{:});
    velnbrcorrCoM = vertcat(velnbrcorrCoM{:});
    accnbrcorr = vertcat(accnbrcorr{:});
    polar_order = vertcat(polar_order{:});
    nematic_order = vertcat(nematic_order{:});
    %% bin distance data
    [~,nbrDistBins,nbrDistbinIdx]  = histcounts(nbrDist,...
        'BinWidth',distBinWidth,'BinLimits',[minDist maxDist]);
    [~,pairDistBins,pairDistbinIdx]  = histcounts(pairdist,...
        'BinWidth',distBinWidth,'BinLimits',[minDist maxDist]);
    % convert bin edges to centres (for plotting)
    nbrDistBins = double(nbrDistBins(1:end-1) + mean(diff(nbrDistBins))/2);
    pairDistBins = double(pairDistBins(1:end-1) + mean(diff(pairDistBins))/2);
    % ignore out of bounds/very small distance values (bin=0) %and bins with only one element
    nbrdistkeepIdcs = nbrDistbinIdx>0;%&ismember(nbrDistbinIdx,find(nbrDistcounts>1));
    % nbrDistBins = nbrDistBins(nbrDistcounts>1);
    nbrDistbinIdx = nbrDistbinIdx(nbrdistkeepIdcs);
    nbrDist = nbrDist(nbrdistkeepIdcs);
    velnbrcorr = velnbrcorr(nbrdistkeepIdcs);
    velnbrcorrSigned = velnbrcorrSigned(nbrdistkeepIdcs);
    velnbrcorrFwd = velnbrcorrFwd(nbrdistkeepIdcs);
    distCoMkeepIdcs = distCoM>minDist;
    distCoM = distCoM(distCoMkeepIdcs);
    velnbrcorrCoM = velnbrcorrCoM(distCoMkeepIdcs);
    accnbrcorr = accnbrcorr(nbrdistkeepIdcs);
    pdistkeepIdcs = pairDistbinIdx>0;%&ismember(pairDistbinIdx,find(pairDistcounts>1));
    % pairDistBins = pairDistBins(pairDistcounts>1);
    dirxcorr = dirxcorr(pdistkeepIdcs);
    velxcorr = velxcorr(pdistkeepIdcs);
    pairDistbinIdx = pairDistbinIdx(pdistkeepIdcs);

    [corr_vn_mean,corr_vn_ci] = grpstats(velnbrcorr,nbrDistbinIdx,{'mean','meanci'});
    [corr_absvn_mean,corr_absvn_ci] = grpstats(abs(velnbrcorr),nbrDistbinIdx,{'mean','meanci'});
    [corr_an_mean,corr_an_ci] = grpstats(accnbrcorr,nbrDistbinIdx,{'mean','meanci'});
    [corr_dir_mean,corr_dir_ci] = grpstats(dirxcorr,pairDistbinIdx,{'mean','meanci'});
    [corr_v_mean,corr_v_ci] = grpstats(velxcorr,pairDistbinIdx,{'mean','meanci'});
    %% plot data
    % correlations
    boundedline(pairDistBins/1000,corr_dir_mean,[corr_dir_mean - corr_dir_ci(:,1), corr_dir_ci(:,2) - corr_dir_mean],...
        'alpha',dircorrFig.Children,'cmap',plotColors(strainCtr,:))
    boundedline(pairDistBins/1000,corr_v_mean,[corr_v_mean - corr_v_ci(:,1), corr_v_ci(:,2) - corr_v_mean],...
        'alpha',velcorrFig.Children,'cmap',plotColors(strainCtr,:))
    boundedline(nbrDistBins/1000,corr_vn_mean,[corr_vn_mean - corr_vn_ci(:,1), corr_vn_ci(:,2) - corr_vn_mean],...
        'alpha',velnbrcorrFig.Children,'cmap',plotColors(strainCtr,:))
    boundedline(nbrDistBins/1000,corr_absvn_mean-1/2,[corr_absvn_mean - corr_absvn_ci(:,1), corr_absvn_ci(:,2) - corr_absvn_mean],...
        'alpha',absvelnbrcorrFig.Children,'cmap',plotColors(strainCtr,:))
    boundedline(nbrDistBins/1000,corr_an_mean,[corr_an_mean - corr_an_ci(:,1), corr_an_ci(:,2) - corr_an_mean],...
        'alpha',accnbrcorrFig.Children,'cmap',plotColors(strainCtr,:))
    
    % 2d histogram of angles between velocity and neighbor, for this
    % strain
    set(0,'CurrentFigure',angleFig)
    anglehist = histogram2conditional(nbrDist/1000,real(acos(velnbrcorr)),2,...
        'XBinLimits',[minDist maxDist]/1000,'BinWidth',[distBinWidth*2/1000 pi/4]); % need to take real part as dot products can be just above or below 1 (less than 1e-6)
    anglehist.EdgeColor = 'none';
    anglehistBinCounts{strainCtr} = anglehist.BinCounts;
    colormap(flipud(cmap_Blues))
    hc = colorbar;
    hc.Label.String = 'relative frequency';
    xlabel('distance to neighbor (mm)')
    ylabel('vel. angle w.r.t. neighbor dir. (rad)')
    ylim([0 pi])
    box on
    figurename = ['figures/correlation/phaseSpecific/anglehist_'  strains{strainCtr} '_' wormnum '_' phase '_data' num2str(dataset) '_' markerType ];     if useJoinedTraj, figurename = [figurename '_jointraj']; end
    exportfig(gcf,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
    
    anglehistSigned = histogram2conditional(nbrDist/1000,real(acos(velnbrcorrSigned)),2,...
        'XBinLimits',[minDist maxDist]/1000,'BinWidth',[distBinWidth*2/1000 pi/4]); % need to take real part as dot products can be just above or below 1 (less than 1e-6)
    anglehistBinCountsSigned{strainCtr} = anglehistSigned.BinCounts;
    
    anglehistFwd = histogram2conditional(nbrDist/1000,real(acos(velnbrcorrFwd)),2,...
        'XBinLimits',[minDist maxDist]/1000,'BinWidth',[distBinWidth*2/1000 pi/4]); % need to take real part as dot products can be just above or below 1 (less than 1e-6)
    anglehistBinCountsFwd{strainCtr} = anglehistFwd.BinCounts;
    
    anglehistCoM = histogram2conditional(distCoM/1000,real(acos(velnbrcorrCoM)),2,...
        'XBinLimits',[minDist 1400+minDist]/1000,'BinWidth',[distBinWidth*2/1000 pi/4]); % need to take real part as dot products can be just above or below 1 (less than 1e-6)
    anglehistBinCountsCoM{strainCtr} = anglehistCoM.BinCounts;
    
    % plot pair-correlation
    pcf = cat(2,pcf{:});
    [lineHandles(strainCtr), ~] = boundedline((distBins(2:end)-distBinWidth/2)/1000,nanmean(pcf,2),...
        [nanstd(pcf,0,2) nanstd(pcf,0,2)]./sqrt(nnz(sum(~isnan(pcf),2))),...
        'alpha',poscorrFig.Children,'cmap',plotColors(strainCtr,:));
    % plot orientational order
    set(0,'CurrentFigure',orderFig)
    subplot(1,2,strainCtr)
    violinplot([polar_order, nematic_order],{'p','n'});
    title(strains{strainCtr})
    ylim([0 1]), box on
    % plot frequency of visited sites
    if  strcmp(wormnum,'40')&& plotDiagnostics
        histogram(visitfreqFig.Children,vertcat(visitfreq{:}),'DisplayStyle','stairs','Normalization','pdf')
    end
end
%% format and export figures
for figHandle = [dircorrFig, velcorrFig, velnbrcorrFig, absvelnbrcorrFig, accnbrcorrFig, poscorrFig] % common formating for all figures
    set(figHandle,'PaperUnits','centimeters')
    figHandle.Children.XLim = [0 maxDist]/1000;
    figHandle.Children.XGrid = 'on';
    figHandle.Children.YGrid = 'on';
    figHandle.Children.Box = 'on';
    set(figHandle.Children,'XTick',distxticks,'XTickLabel',num2str(distxticks'))
    if figHandle~=poscorrFig&&figHandle~=absvelnbrcorrFig
        figHandle.Children.YLim = [-0.5 0.5];
    end
end

% directional correlation
ylabel(dircorrFig.Children,'orientational correlation')
xlabel(dircorrFig.Children,'distance between pair (mm)')
legend(dircorrFig.Children,lineHandles,strains)
figurename = ['figures/correlation/phaseSpecific/dircrosscorr_' wormnum '_' phase '_data' num2str(dataset) '_' markerType ];     if useJoinedTraj, figurename = [figurename '_jointraj']; end
exportfig(dircorrFig,[figurename '.eps'],exportOptions)
system(['epstopdf ' figurename '.eps']);
system(['rm ' figurename '.eps']);

% velocity correlation
ylabel(velcorrFig.Children,'velocity correlation')
xlabel(velcorrFig.Children,'distance between pair (mm)')
legend(velcorrFig.Children,lineHandles,strains)
figurename = ['figures/correlation/phaseSpecific/velcrosscorr_' wormnum '_' phase '_data' num2str(dataset) '_' markerType ];     if useJoinedTraj, figurename = [figurename '_jointraj']; end
exportfig(velcorrFig,[figurename '.eps'],exportOptions)
system(['epstopdf ' figurename '.eps']);
system(['rm ' figurename '.eps']);

% correlation of velocity with direction to neighbor
ylabel(velnbrcorrFig.Children,'vel. to nbr dir. corr.')
xlabel(velnbrcorrFig.Children,'distance to neighbor (mm)')
legend(velnbrcorrFig.Children,lineHandles,strains)
figurename = ['figures/correlation/phaseSpecific/velnbrcorr_' wormnum '_' phase '_data' num2str(dataset) '_' markerType ];     if useJoinedTraj, figurename = [figurename '_jointraj']; end
exportfig(velnbrcorrFig,[figurename '.eps'],exportOptions)
system(['epstopdf ' figurename '.eps']);
system(['rm ' figurename '.eps']);

% absolute value of correlation of velocity with direction to neighbor
ylabel(absvelnbrcorrFig.Children,'$\langle|\vec{v}\cdot\vec{r}_\mathrm{nbr}|/(|\vec{v}||\vec{r}_\mathrm{nbr}|)\rangle - \frac{1}{2}$','Interpreter','LaTeX')
xlabel(absvelnbrcorrFig.Children,'distance to neighbor (mm)')
absvelnbrcorrFig.Children.YLim = [0 1];
legend(absvelnbrcorrFig.Children,lineHandles,strains)
figurename = ['figures/correlation/phaseSpecific/absvelnbrcorr_' wormnum '_' phase '_data' num2str(dataset) '_' markerType ];     if useJoinedTraj, figurename = [figurename '_jointraj']; end
exportfig(absvelnbrcorrFig,[figurename '.eps'],exportOptions)
system(['epstopdf ' figurename '.eps']);
system(['rm ' figurename '.eps']);

% correlation of acceleration with direction to neighbor
ylabel(accnbrcorrFig.Children,'acc.-dir. to neighbor corr.')
xlabel(accnbrcorrFig.Children,'distance to neighbor (mm)')
legend(accnbrcorrFig.Children,lineHandles,strains)
figurename = ['figures/correlation/phaseSpecific/accnbrcorr_' wormnum '_' phase '_data' num2str(dataset) '_' markerType ];     if useJoinedTraj, figurename = [figurename '_jointraj']; end
exportfig(accnbrcorrFig,[figurename '.eps'],exportOptions)
system(['epstopdf ' figurename '.eps']);
system(['rm ' figurename '.eps']);

% 2d histogram of angles between velocity and neighbor, difference between strains
load('~/Dropbox/Utilities/colormaps_ascii/diverging/cmap_BuRd.txt')
for plotCtr =1:4
    switch plotCtr
        case 1
            binCountsToPlot = anglehistBinCounts;
            histhandle = anglehist;
            thisYLabel = 'vel. angle w.r.t. nbr dir. (rad)';
            figurename = ['figures/correlation/phaseSpecific/anglehistDiff_'  strains{strainCtr} '_' wormnum '_' phase '_data' num2str(dataset) '_' markerType ];     if useJoinedTraj, figurename = [figurename '_jointraj']; end
        case 2
            binCountsToPlot = anglehistBinCountsSigned;
            histhandle = anglehist;
            thisYLabel = 'signed vel. angle w.r.t. nbr dir. (rad)';
            figurename = ['figures/correlation/phaseSpecific/anglehistDiffSigned_'  strains{strainCtr} '_' wormnum '_' phase '_data' num2str(dataset) '_' markerType ];     if useJoinedTraj, figurename = [figurename '_jointraj']; end
        case 3
            binCountsToPlot = anglehistBinCountsFwd;
            histhandle = anglehist;
            thisYLabel = 'fwd vel. angle w.r.t. nbr dir. (rad)';
            figurename = ['figures/correlation/phaseSpecific/anglehistDiffFwd_'  strains{strainCtr} '_' wormnum '_' phase '_data' num2str(dataset) '_' markerType ];     if useJoinedTraj, figurename = [figurename '_jointraj']; end
        case 4
            binCountsToPlot = anglehistBinCountsCoM;
            histhandle = anglehistCoM;
            thisYLabel = 'vel. angle w.r.t. nbr CoM (rad)';
            figurename = ['figures/correlation/phaseSpecific/anglehistDiffCoM_'  strains{strainCtr} '_' wormnum '_' phase '_data' num2str(dataset) '_' markerType ];     if useJoinedTraj, figurename = [figurename '_jointraj']; end
    end
    angleDiffFig = figure;
    set(0,'CurrentFigure',angleDiffFig)
    % to plot histogram with negative values, offset bincounts
    bincountOffset = max(binCountsToPlot{2}(:));
    histogram2('XBinEdges',histhandle.XBinEdges,'YBinEdges',histhandle.YBinEdges,...
        'BinCounts',(binCountsToPlot{1} - binCountsToPlot{2} + bincountOffset),'DisplayStyle','tile','EdgeColor','none')
    colormap(flipud(cmap_BuRd))
    hc = colorbar;
    caxisRaw = caxis;
    caxis(bincountOffset + [-1 1]*max(abs(caxisRaw - bincountOffset))) % center color range around zero
    caxisCentered = caxis;
    hc.Ticks = linspace(caxisCentered(1),caxisCentered(2),7);
    hc.TickLabels = num2str(hc.Ticks' - bincountOffset,2); % correct labels for offset
    hc.Label.String = 'difference in rel. freq., npr-1 to N2';
    xlabel('distance to neighbor (mm)')
    ylim([0 pi])
%     xlim([0 maxDist/1000])
    box on
    ylabel(thisYLabel)
    exportfig(gcf,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end

% pair-correlation function
poscorrFig.Children.YLim(1) = 0;
% poscorrFig.Children.YTick = 0:2:round(poscorrFig.Children.YLim(2));
ylabel(poscorrFig.Children,'pair correlation')
xlabel(poscorrFig.Children,'distance r (mm)')
legend(poscorrFig.Children,lineHandles,strains)
figurename = ['figures/correlation/phaseSpecific/paircorrelationfunction_' wormnum '_' phase '_data' num2str(dataset) '_' markerType];
if useJoinedTraj, figurename = [figurename '_jointraj']; end
exportfig(poscorrFig,[figurename '.eps'],exportOptions)
system(['epstopdf ' figurename '.eps']);
system(['rm ' figurename '.eps']);

% polar and nematic order
figurename = ['figures/correlation/phaseSpecific/polarNematicOrder_' wormnum '_' phase '_data' num2str(dataset) '_' markerType];
if useJoinedTraj, figurename = [figurename '_jointraj']; end
exportfig(orderFig,[figurename '.eps'],exportOptions)
system(['epstopdf ' figurename '.eps']);
system(['rm ' figurename '.eps']);
% heatmap of sites visited
if  strcmp(wormnum,'40')&& plotDiagnostics
    visitfreqFig.Children.XScale = 'log';
    visitfreqFig.Children.YScale = 'log';
    %         visitfreqFig.Children.XLim = [4e-5 1e-1];
    xlabel(visitfreqFig.Children,'site visit frequency, f')
    ylabel(visitfreqFig.Children,'pdf p(f)')
    legend(visitfreqFig.Children,strains)
    figurename = ['figures/visitfreq_' wormnum '_' markerType];
    exportfig(visitfreqFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end
