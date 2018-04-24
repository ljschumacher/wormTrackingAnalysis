% dataset: 1 or 2. To specify which dataset to run the script for.
% phase: 'joining', 'fullMovie', or 'sweeping'. Script defines stationary phase as: starts at 10% into the movie, and stops at 60% into the movie (HA and N2) or at specified stopping frames (npr-1).
% wormnum: '40', or 'HD'
% markerType: 'pharynx', or 'bodywall'
% plotDiagnostics: true (default) or false

close all
clear

dataset = 2;
phase = 'sweeping';
wormnum = '40';
markerType = 'pharynx';
plotDiagnostics = true;

addpath('auxiliary/')
addpath('visualisation/')

%% set fixed parameters

if dataset ==1
    strains = {'npr1','N2'};%{'npr1','HA','N2'}
    assert(~strcmp(markerType,'bodywall'),'Bodywall marker for dataset 1 not available')
elseif dataset ==2
    strains = {'npr1','N2'};
end
useJoinedTraj = true;

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
% plotting parameters
plotColors = lines(nStrains);
if plotDiagnostics, visitfreqFig = figure; hold on, end
% export fig parameters
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

%% loop through strains
for strainCtr = 1:nStrains
    %% load file lists
    if dataset == 1
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_list.xlsx'],1,'A1:E15','basic');
    elseif dataset == 2
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_' channelStr '_list.xlsx'],1,'A1:E15','basic');
    end
    if ~useJoinedTraj
        filenames = strrep(filenames,'/data2/shared/data/twoColour/Results/',...
            '/end/home/lschumac/databackup/data/twoColour/ResultsUnjoinedTrajectories/');
    end
    numFiles = length(filenames);
    if strcmp(wormnum,'40'), visitfreq = cell(numFiles,1); end
    %% loop through files
    for fileCtr = 1:numFiles % can be parfor
        filename = filenames{fileCtr};
        %% load tracking data
        trajData = h5read(filename,'/trajectories_data');
        blobFeats = h5read(filename,'/blob_features');
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

        %% heat map of sites visited - this only makes sense for 40 worm dataset where we don't move the camera
        if strcmp(wormnum,'40')&& plotDiagnostics
            siteVisitFig = figure;
            x = trajData.coord_x(trajData.filtered);
            y = trajData.coord_y(trajData.filtered);
            h=histogram2(x*pixelsize/1000,y*pixelsize/1000,...
                'DisplayStyle','tile','EdgeColor','none','Normalization','pdf');
            visitfreq{fileCtr} = h.Values(:);
            cb = colorbar; cb.Label.String = '# visited';
            axis equal
            xlabel('x (mm)'), ylabel('y (mm)')
            title([strains{strainCtr} ' ' strrep(filename(end-32:end-18),'/','')])
            set(siteVisitFig,'PaperUnits','centimeters')
            figurename = ['figures/individualRecordings/' strains{strainCtr}...
                '_' strrep(strrep(filename(end-32:end-18),' ',''),'/','') '_sitesVisited' '_' phase '_data' num2str(dataset)];
            %exportfig(siteVisitFig,[figurename '.eps'],exportOptions)
        end
    end
     %% combine data from multiple files
         % plot frequency of visited sites
    if  strcmp(wormnum,'40')&& plotDiagnostics
        histogram(visitfreqFig.Children,vertcat(visitfreq{:}),'DisplayStyle','stairs','Normalization','pdf')
    end
end
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
end