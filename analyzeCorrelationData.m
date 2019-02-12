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
maxSpeed = 1000;
maxDist = 2000;
distBins = 0:distBinWidth:maxDist;
% define functions for grpstats
% mad1 = @(x) mad(x,1); % median absolute deviation
% % alternatively could use boxplot-style confidence intervals on the mean,
% % which are 1.57*iqr/sqrt(n) - unclear how justified this is
% iqrci = @(x) 1.57*iqr(x)/sqrt(numel(x));
% % or one could use a bootstrapped confidence interval
bootserr = @(x) bootci(1e2,{@median,x},'alpha',0.05,'Options',struct('UseParallel',false));
% plotting parameters
plotColors = lines(nStrains);
if plotDiagnostics, visitfreqFig = figure; hold on, end
dircorrxticks = 0:0.5:(maxDist/1000);
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

poscorrFig = figure; hold on
orderFig = figure; hold on
lineHandles = NaN(nStrains,1);
%% loop through strains
for strainCtr = 1:nStrains
    dircorrFig = figure; hold on
    velcorrFig = figure; hold on
    velnbrcorrFig = figure; hold on
    accnbrcorrFig = figure; hold on
    %% load file lists
    if dataset == 1
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_list_lslx.xlsx'],1,'A1:E15','basic');
    elseif dataset == 2
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_' channelStr '_list_lslx.xlsx'],1,'A1:E15','basic');
    end
    if ~useJoinedTraj
        filenames = strrep(filenames,'/data2/shared/data/twoColour/Results/',...
            '/end/home/lschumac/databackup/data/twoColour/ResultsUnjoinedTrajectories/');
    end
    numFiles = length(filenames);
    if strcmp(wormnum,'40'), visitfreq = cell(numFiles,1); end
    %% intialise variables to store results
    dirxcorr = cell(numFiles,1); % for calculating directional cross-correlation
    velxcorr = cell(numFiles,1); % for calculating velocity cross-correlation
    velnbrcorr = cell(numFiles,1); % for calculating velocity correlation with nearest neighbour position
    accnbrcorr = cell(numFiles,1);
    pairdist = cell(numFiles,1);
    nbrDist= cell(numFiles,1); % for the weighted neighbour distances
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
        assert(size(skelData,2)==2,['Wrong skeleton size for ' filename])
        assert(size(skelData,3)==length(trajData.frame_number),['Wrong number of skeleton frames for ' filename])
        assert(length(blobFeats.velocity_x)==length(trajData.frame_number)&&...
            length(blobFeats.signed_speed)==length(trajData.frame_number),['Wrong number of speed frames for ' filename])
        if all(isnan(skelData(:))), warning(['all skeleton are NaN for ' filename]),end
        %% randomly sample which frames to analyze
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr);
        numFrames = round((lastFrame-firstFrame)/frameRate/3);
        framesAnalyzed = randperm((lastFrame-firstFrame),numFrames) + firstFrame; % randomly sample frames without replacement
        %% filter data for worms
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
        accnbrcorr{fileCtr} = cell(numFrames,1);
        pairdist{fileCtr} = cell(numFrames,1);
        nbrDist{fileCtr}= cell(numFrames,1);
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
                % assume speeds>maxSpeed are tracking errors
                speedOutlierLogIdcs = sqrt(vx.^2 + vy.^2)>maxSpeed;
                if any(speedOutlierLogIdcs)
                    vx(speedOutlierLogIdcs) = NaN;
                    vy(speedOutlierLogIdcs) = NaN;
                    accx(speedOutlierLogIdcs) = NaN;
                    accy(speedOutlierLogIdcs) = NaN;
                end
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
                %% calculate direction towards neighbour(s) for each worm
                dx = (x - x')*pixelsize; dy = (y - y')*pixelsize;
                dD = sqrt(dx.^2 + dy.^2);
                % keep only nbr-distances below maxDist
                nbrKeepLogIdcs = dD(:)<=maxDist;
                
                nbrDist{fileCtr}{frameCtr} = dD(nbrKeepLogIdcs);
                velnbrcorr{fileCtr}{frameCtr} = reshape(double(vectorPairedCorrelation2D(vx,vy,dx./dD,dy./dD,true,false)),[],1); % reshape to column-vector
                velnbrcorr{fileCtr}{frameCtr} = velnbrcorr{fileCtr}{frameCtr}(nbrKeepLogIdcs);
                % calculate correlation between change in direction and nbr positions
                accnbrcorr{fileCtr}{frameCtr} = reshape(double(vectorPairedCorrelation2D(accx,accy,dx./dD,dy./dD,true,false)),[],1); % reshape to column-vector
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
        accnbrcorr{fileCtr} = vertcat(accnbrcorr{fileCtr}{:});
        pairdist{fileCtr} = horzcat(pairdist{fileCtr}{:});
        nbrDist{fileCtr} = vertcat(nbrDist{fileCtr}{:});
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
    pairdist = horzcat(pairdist{:});
    dirxcorr = horzcat(dirxcorr{:});
    velxcorr = horzcat(velxcorr{:});
    velnbrcorr = vertcat(velnbrcorr{:});
    accnbrcorr = vertcat(accnbrcorr{:});
    polar_order = vertcat(polar_order{:});
    nematic_order = vertcat(nematic_order{:});
    
    %% plot data
    % correlations
    set(0,'CurrentFigure',dircorrFig)
    histogram2conditional(pairdist/1000,dirxcorr,2)
    set(0,'CurrentFigure',velcorrFig)
    histogram2conditional(pairdist/1000,velxcorr,2)
    set(0,'CurrentFigure',velnbrcorrFig)
    histogram2conditional(nbrDist/1000,velnbrcorr,2)
    set(0,'CurrentFigure',accnbrcorrFig)
    histogram2conditional(nbrDist/1000,accnbrcorr,2)
    
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
    %% format and export figures for this strain
    for figHandle = [dircorrFig, velcorrFig, velnbrcorrFig, accnbrcorrFig] % common formating for all figures
        set(figHandle,'PaperUnits','centimeters')
    end
    % directional correlation
    dircorrFig.Children.YLim = [-1 1];
    dircorrFig.Children.XLim = [0 maxDist]/1000;
    set(dircorrFig.Children,'XTick',dircorrxticks,'XTickLabel',num2str(dircorrxticks'))
    ylabel(dircorrFig.Children,'orientational correlation')
    xlabel(dircorrFig.Children,'distance between pair (mm)')
    figurename = ['figures/correlation/phaseSpecific/dircrosscorr_' strains{strainCtr} '_' wormnum '_' phase '_data' num2str(dataset) '_' markerType ];     if useJoinedTraj, figurename = [figurename '_jointraj']; end
    exportfig(dircorrFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
    % velocity correlation
    velcorrFig.Children.YLim = [-1 1];
    velcorrFig.Children.XLim = [0 maxDist]/1000;
    set(velcorrFig.Children,'XTick',dircorrxticks,'XTickLabel',num2str(dircorrxticks'))
    ylabel(velcorrFig.Children,'velocity correlation')
    xlabel(velcorrFig.Children,'distance between pair (mm)')
    figurename = ['figures/correlation/phaseSpecific/velcrosscorr_' strains{strainCtr} '_' wormnum '_' phase '_data' num2str(dataset) '_' markerType ];     if useJoinedTraj, figurename = [figurename '_jointraj']; end
    exportfig(velcorrFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
    % correlation of velocity with direction to neighbour
    velnbrcorrFig.Children.YLim = [-1 1];
    velnbrcorrFig.Children.XLim = [0 maxDist]/1000;
    set(velnbrcorrFig.Children,'XTick',dircorrxticks,'XTickLabel',num2str(dircorrxticks'))
    ylabel(velnbrcorrFig.Children,'velocity-nbr corr.')
    xlabel(velnbrcorrFig.Children,'neighbour distance (mm)')
    figurename = ['figures/correlation/phaseSpecific/velnbrcorr_' strains{strainCtr} '_' wormnum '_' phase '_data' num2str(dataset) '_' markerType ];     if useJoinedTraj, figurename = [figurename '_jointraj']; end
    exportfig(velnbrcorrFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
    % correlation of acceleration with direction to neighbour
    accnbrcorrFig.Children.YLim = [-1 1];
    accnbrcorrFig.Children.XLim = [0 maxDist]/1000;
    accnbrcorrFig.Children.Box = 'on';
    set(accnbrcorrFig.Children,'XTick',dircorrxticks,'XTickLabel',num2str(dircorrxticks'))
    ylabel(accnbrcorrFig.Children,'acceleration-nbr corr.')
    xlabel(accnbrcorrFig.Children,'neighbour distance (mm)')
    figurename = ['figures/correlation/phaseSpecific/accnbrcorr_' strains{strainCtr} '_' wormnum '_' phase '_data' num2str(dataset) '_' markerType ];     if useJoinedTraj, figurename = [figurename '_jointraj']; end
    exportfig(accnbrcorrFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end
%% format and export figures
for figHandle = [poscorrFig] % common formating for all figures
    set(figHandle,'PaperUnits','centimeters')
end
% pair-correlation function
poscorrFig.Children.YLim(1) = 0;
poscorrFig.Children.XLim = [0 maxDist]/1000;
poscorrFig.Children.XTick = 0:0.5:(maxDist/1000);
% poscorrFig.Children.YTick = 0:2:round(poscorrFig.Children.YLim(2));
poscorrFig.Children.Box = 'on';
poscorrFig.Children.XGrid = 'on';
poscorrFig.Children.YGrid = 'on';
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
