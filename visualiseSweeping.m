% dataset: 1 or 2. To specify which dataset to run the script for.
% phase: 'joining', 'fullMovie', or 'sweeping'. Script defines stationary phase as: starts at 10% into the movie, and stops at 60% into the movie (HA and N2) or at specified stopping frames (npr-1).
% wormnum: '40', or 'HD'
% markerType: 'pharynx', or 'bodywall'
% plotDiagnostics: true (default) or false

close all
clear

dataset = 2;
phase = 'fullMovie';
wormnum = '40';
markerType = 'pharynx';
numMovieSlices = 10;

if numMovieSlices>1
    % frame slice: n slices of 1 minutes (9 fps * 60s)
    frameSlices = zeros(numMovieSlices,2);
    for sliceCtr = 1:numMovieSlices
        frameSlices(sliceCtr,:) = [round(32400/numMovieSlices*(sliceCtr-1)) round(32400/numMovieSlices*(sliceCtr-1)+540)];
    end
end


addpath('auxiliary/')
addpath('visualisation/')

%% set fixed parameters

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
else
    error('unknown marker type specified, should be pharynx or bodywall')
end
% analysis parameters
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
% export fig parameters
exportOptions = struct('Format','EPS2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

%% loop through strains
for strainCtr = 1%:nStrains
    %% load file lists
    if dataset == 1
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_list.xlsx'],1,'A1:E15','basic');
    elseif dataset == 2
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_' channelStr '_list.xlsx'],1,'A1:E15','basic');
    end
    numFiles = length(filenames);
    if strcmp(wormnum,'40'), visitfreq = cell(numFiles,1); end
    %% loop through files
    for fileCtr = 1%:numFiles % can be parfor
        filename = filenames{fileCtr};
        %% make new video
        writerObj = VideoWriter(['figures/individualRecordings/' strains{strainCtr} '_' strrep(strrep(filename(end-32:end-18),' ',''),'/','') '_sitesVisited.avi']);
        writerObj.FrameRate = 2;
        open(writerObj);
        %% load tracking data
        trajData = h5read(filename,'/trajectories_data');
        blobFeats = h5read(filename,'/blob_features');
        %% sample which frames to analyze
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr);
        %% filter data for worms
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
        % select frame slices
        for sliceCtr = 1:size(frameSlices,1)
            sliceStart = frameSlices(sliceCtr,1);
            sliceEnd = frameSlices(sliceCtr,2);
            sliceLogInd = false(size(trajData.filtered));
            for frameCtr = sliceStart:sliceEnd
                sliceLogInd(frameCtr==trajData.frame_number)=true;
            end
            %% heat map of sites visited - this only makes sense for 40 worm dataset where we don't move the camera
            if strcmp(wormnum,'40')
                siteVisitFig = figure;
                x = trajData.coord_x(trajData.filtered & sliceLogInd);
                y = trajData.coord_y(trajData.filtered & sliceLogInd);
                h=histogram2(x*pixelsize/1000,y*pixelsize/1000,...
                    'DisplayStyle','tile','EdgeColor','none','Normalization','count');
                visitfreq{fileCtr} = h.Values(:);
                cb = colorbar; cb.Label.String = '# visited';
                caxis([0 1500]);
                xlim([0 12]);
                ylim([0 12]);
                xlabel('x (mm)'), ylabel('y (mm)')
                title([strains{strainCtr} ' ' strrep(filename(end-32:end-18),'/','') ', ' num2str(round((sliceCtr-1)*60/numMovieSlices)) 'min'])
                set(siteVisitFig,'PaperUnits','centimeters')
                figurename = ['figures/individualRecordings/' strains{strainCtr}...
                    '_' strrep(strrep(filename(end-32:end-18),' ',''),'/','') '_sitesVisited' '_' phase '_slice' num2str(sliceCtr) '_data' num2str(dataset)];
                saveas(gcf,[figurename '.tif'])
                % write heatmap to video
                image = imread([figurename '.tif']);
                frame = im2frame(image);
                writeVideo(writerObj,frame)
                system(['rm ' figurename '.tif']);
            end
        end
        close(writerObj);
    end
end