function [] = analyzeAutoCorrelation(dataset,phase,wormnum,markerType,plotDiagnostics)
% calculate velocity auto correlation, to parameterise noise
% INPUTS
% dataset: 1 or 2. To specify which dataset to run the script for.
% phase: 'joining', 'fullMovie', or 'sweeping'. Script defines stationary phase as: starts at 10% into the movie, and stops at 60% into the movie (HA and N2) or at specified stopping frames (npr-1).
% wormnum: '40', or 'HD'
% markerType: 'pharynx', or 'bodywall'
% OUTPUTS
% none returned, but figures are exported
% issues/to-do:
% - correct for uneven spacing of frames

addpath('auxiliary/')
addpath('filters/')

%% set fixed parameters

if nargin<5
    plotDiagnostics = false; % true or false
    if nargin<4
        markerType = 'bodywall';
    end
end

if dataset ==1
    strains = {'npr1','N2'};%{'npr1','HA','N2'}
    assert(~strcmp(markerType,'bodywall'),'Bodywall marker for dataset 1 not available')
elseif dataset ==2
    strains = {'npr1','N2'};
end

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
maxSpeed = 500;
minSpeedPerFrame = pixelsize/2; % assume these worms to non-moving, and exclude from velocity correlation analysis (as that will be inaccurate for very low magnitude velocities)
maxLag = 30; % max lag in seconds

% plotting parameters
plotColors = lines(nStrains);

% export fig parameters
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);
%% initialize figures
vaccorrFig = figure; hold on

%% temporarily skip N2 for 40 worm bodywall data as those trajectories have not been curated
if strcmp(wormnum,'40')&&strcmp(markerType,'bodywall')
    nStrains = 1;
end
%% loop through strains
for strainCtr = 1:nStrains
    %% load file lists
    if dataset == 1
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_list_lslx.xlsx'],1,'A1:E15','basic');
    elseif dataset == 2
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_' channelStr '_list_lslx.xlsx'],1,'A1:E15','basic');
    end
    
    numFiles = length(filenames);
    %% intialise variables to store results
    velacorr = cell(numFiles,1); % for calculating velocity auto-correlation
    wormFrameCounts = cell(numFiles,1);
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
        %% select subset of frames to analyze
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr);
        framesAnalyzed = firstFrame:lastFrame;
        numFrames = nnz(framesAnalyzed);
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
            if strcmp(wormnum,'1W')
                trajData.filtered = trajData.filtered&...
                    filterSkelLength(skelData,pixelsize,minSkelLength,maxSkelLength);
            else
                trajData.filtered = trajData.filtered&logical(trajData.is_good_skel)&...
                    filterSkelLength(skelData,pixelsize,minSkelLength,maxSkelLength);
            end
        end
        % apply phase restriction
        phaseFilter_logInd = trajData.frame_number < lastFrame & trajData.frame_number > firstFrame;
        trajData.filtered(~phaseFilter_logInd)=false;
        %% calculate stats
        if ~isfield(trajData,'worm_index_manual')
            trajData.worm_index_manual = trajData.worm_index_joined;
            warning(['working with non-curated worm indeces for ' filename])
        end
        wormIDs = unique(trajData.worm_index_manual(trajData.filtered));
        numWorms = numel(wormIDs);
        % check how many frames we have or each worm
        for wormCtr = 1:numWorms
            wormID = wormIDs(wormCtr);
            thisWormLogInd = trajData.worm_index_manual==wormID&trajData.filtered;
            % if the worm doesn't have enough frames, don't consider it
            if nnz(thisWormLogInd)<maxLag*frameRate
                trajData.filtered(thisWormLogInd) = false;
            end
        end
        % reset worm IDs after filtering for trajectoru length
        wormIDs = unique(trajData.worm_index_manual(trajData.filtered));
        numWorms = numel(wormIDs);
        if numWorms>0
            % initialise variables to store results for this file
            velacorr{fileCtr} = NaN(numWorms,1+maxLag*frameRate); % for calculating velocity cross-correlation
            wormFrameCounts{fileCtr} = NaN(numWorms,1);
            % calculate speeds
            if strcmp(markerType,'bodywall')
                [ v, velocities_x, velocities_y, ~, ~, ~ ] = ...
                    calculateSpeedsFromSkeleton(trajData,skelData,1:5,...
                    pixelsize,frameRate,true,0); % last arg is smoothing window
            end
            %          filter speed outliers
            %                         speedOutlierLogIdcs = v>maxSpeed|v<minSpeed;
            %         velocities_x(speedOutlierLogIdcs) = NaN;
            %         velocities_y(speedOutlierLogIdcs) = NaN;
            
            for wormCtr = 1:numWorms
                wormID = wormIDs(wormCtr);
                thisWormLogInd = trajData.worm_index_manual==wormID&trajData.filtered;
                %             % shift velocity time-series by frame-number (so that
                %             % deltaFrame between consecutive entries is 1 for the vac calc)
                %             thisWormFrameNumbers = trajData.frame_number(thisWormLogInd);
                %             thisWormFirstFrame = min(thisWormFrameNumbers);
                %             thisWormLastFrame = max(thisWormFrameNumbers);
                %             thisWormVelShifted = NaN(2,thisWormLastFrame - thisWormFirstFrame + 1);
                %             for frame = thisWormFrameNumbers'
                %                 thisWormVelShifted(:,frame - thisWormFirstFrame + 1) = ...
                %                     [velocities_x(thisWormLogInd&frame==trajData.frame_number);
                %                     velocities_y(thisWormLogInd&frame==trajData.frame_number)]
                %             end
                % check frames are contiguous
                % %             assert(1==unique(diff(trajData.frame_number(thisWormLogInd))),'worm frames are not contiguous')
                % calculate autocorrelation
                fullVac = vectorAutoCorrelation(...
                    [velocities_x(thisWormLogInd)'; velocities_y(thisWormLogInd)'],...
                    maxLag*frameRate,true,true);
                for lagCtr = 1:(maxLag*frameRate+1)
                    velacorr{fileCtr}(wormCtr,lagCtr) = nanmean(fullVac{lagCtr});
                end
                % keep track of how many frames we have on this worm
                wormFrameCounts{fileCtr}(wormCtr) = nnz(thisWormLogInd);
            end
        else
            warning(['no worms with sufficiently long trajectories for '...
                filename])
        end
    end
    %% combine data from multiple files
    velacorr = vertcat(velacorr{:});
    wormFrameCounts = vertcat(wormFrameCounts{:});
    
    %% plot data
    plot(vaccorrFig.Children,linspace(0,maxLag,maxLag*frameRate+1),...
        weightedStats(velacorr,wormFrameCounts,'w'),...
        'Color',plotColors(strainCtr,:),'LineWidth',2)
end
%% format and export figures
for figHandle = [vaccorrFig] % common formating for all figures
    set(figHandle,'PaperUnits','centimeters')
    figHandle.Children.XLim = [0 maxLag];
    figHandle.Children.XGrid = 'on';
    figHandle.Children.YGrid = 'on';
    figHandle.Children.Box = 'on';
    figHandle.Children.YLim = [0 1];
end

% velocity correlation
ylabel(vaccorrFig.Children,'velocity auto-correlation')
xlabel(vaccorrFig.Children,'lag (s)')
legend(vaccorrFig.Children,strains)
figurename = ['figures/correlation/phaseSpecific/velautocorr_' wormnum '_' phase '_data' num2str(dataset) '_' markerType ]; 
exportfig(vaccorrFig,[figurename '.eps'],exportOptions)
system(['epstopdf ' figurename '.eps']);
system(['rm ' figurename '.eps']);
end
