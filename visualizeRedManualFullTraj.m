
clear
close all

%% set parameters
phase = 'joining'; % 'fullMovie', 'joining', or 'sweeping'.
dataset = 2; % 1 or 2
marker = 'bodywall'; % 'pharynx' or 'bodywall'
strains = {'npr1'}; % {'npr1','N2'}
wormnums = {'40'};% {'40'};
wormcats = {'leaveCluster','loneWorm'}; %'leaveCluster','loneWorm'
smoothing = false;
postExitDuration = 5; % duration (in seconds) after a worm exits a cluster to be included in the leave cluster analysis
minTrajDuration = 1; % duration (in seconds) of minimum traj length
maxTrajDuration = 5;  % duration (in seconds) of maximum traj length % may set to 1.5 to truncate loneWorm traj to match those of leaveCluster traj length
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
saveResults = true;
visualiseSampleTraj = true; % true or false

if visualiseSampleTraj == true
    featureToSample = 'headAngTotal'; % 'headAngTotal','headAngNorm', or 'headAngSpeed'
    if strcmp(featureToSample,'headAngTotal') | strcmp(featureToSample,'headAngSpeed')
        headAngRanges = [0, 0.25; pi/2-0.25, pi/2+0.25; pi-0.25, pi+0.25; 3/2*pi-0.25, 3/2*pi+0.25; 2*pi-0.25, 2*pi];
    elseif strcmp(featureToSample,'headAngNorm')
        headAngRanges = [0 0.01; 0.01 0.02; 0.02 0.03; 0.03 0.05; 0.05 1];
    else
        warning('Wrong feature selected for trajectory visualisation')
    end
end

maxBlobSize_r = 2.5e5;
minSkelLength_r = 850;
maxSkelLength_r = 1500;
minNeighbrDist = 2000;
inClusterNeighbourNum = 3;

load('exportOptions.mat')

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        % load file list
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_r_list.xlsx'],1,'A1:E15','basic');
        numFiles = length(filenames);
        % create empty cell arrays to hold individual file values, so they can be pooled for a given strain/density combination
        for wormcatCtr = 1:length(wormcats)
            headAngTotal.(wormcats{wormcatCtr}) = cell(numFiles,1);
            headAngNorm.(wormcats{wormcatCtr}) = cell(numFiles,1);
            headAngSpeed.(wormcats{wormcatCtr}) = cell(numFiles,1);
            frameRunLengths.(wormcats{wormcatCtr}) = cell(numFiles,1);
        end
        
        if visualiseSampleTraj
            % save up to 500 sets of xy coordinates for paths that fall within a certain angular speed range to plot sample trajectories
            for wormcatCtr = 1:length(wormcats)
                headAngSampleTraj.(wormcats{wormcatCtr}) = cell(500,2,size(headAngRanges,1));
                headAngSampleCtr.(wormcats{wormcatCtr}) = ones(size(headAngRanges,1),1);
            end
        end
        
        %% go through individual movies
        for fileCtr = 1:numFiles
            %% load data
            filename = filenames{fileCtr}
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton'); % in pixels
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            features = h5read(strrep(filename,'skeletons','feat_manual'),'/features_timeseries');
            
            %% obtain all xy coordinates
            worm_xcoords = squeeze(skelData(1,:,:))'*pixelsize; % turn pixel to microns
            worm_ycoords = squeeze(skelData(2,:,:))'*pixelsize;
             % restrict to head nodes only (8 out of 49 nodes) for bodywall data
            worm_xcoords = worm_xcoords(:,1:8);
            worm_ycoords = worm_ycoords(:,1:8);
                   
            %% filter worms by various criteria
            % filter red by manually joined traj
            trajData.filtered = ismember(trajData.worm_index_manual,int32(features.worm_index));
            % filter red by blob size and intensity
            if contains(filename,'55')||contains(filename,'54')
                intensityThreshold = 80;
            else
                intensityThreshold = 40;
            end
            trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                intensityThreshold,maxBlobSize_r);
            % filter red by skeleton length
            trajData.filtered = trajData.filtered&logical(trajData.is_good_skel)...
                &filterSkelLength(skelData,pixelsize,minSkelLength_r,maxSkelLength_r);
            % apply phase restriction
            [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr);
            phaseFrameLogInd = trajData.frame_number <= lastFrame & trajData.frame_number >= firstFrame;
            trajData.filtered(~phaseFrameLogInd) = false;
            % find worms that have just left a cluster vs lone worms
            [leaveClusterLogInd, loneWormLogInd,~,~] = findWormCategory(filename,inClusterNeighbourNum,minNeighbrDist,postExitDuration);
            
            %% plot worm full vs. categorised traj's
            % go worm by worm and plot full traj for the specified phase
            uniqueWorms = unique(features.worm_index);
            for wormCtr = 1:numel(uniqueWorms)
                uniqueWorm = uniqueWorms(wormCtr);
                wormLogInd = trajData.worm_index_manual==uniqueWorms(wormCtr) ;
                wormtraj_xcoords = worm_xcoords(wormLogInd);
                wormtraj_ycoords = worm_ycoords(wormLogInd);
                fullTraj = figure; hold on 
                plot(wormtraj_xcoords,wormtraj_ycoords,'b')
                
                % now plot individually categorised traj on top of the full traj
                wormLogInd = trajData.worm_index_manual==uniqueWorms(wormCtr) & trajData.filtered;
                for wormcatCtr = 1:length(wormcats)
                    wormCatLogInd = wormLogInd & eval([wormcats{wormcatCtr} 'LogInd']);
                    wormFrames = trajData.frame_number(wormCatLogInd)';
                    
                    % break down frames into continuous trajectories
                    if ~isempty(wormFrames)
                        continuousFrameRuns = getContinuousTraj(wormFrames);
                        
                        % filter for minimum traj length
                        continuousFramesMinLengthLogInd = cellfun(@numel,continuousFrameRuns)>=round(minTrajDuration*frameRate);
                        continuousFrameRuns = continuousFrameRuns(continuousFramesMinLengthLogInd);
                        
                        % go through each traj that fits the min length criteria to obtain xy coordinates
                        if size(continuousFrameRuns,2)>100
                            warning('more trajectories present than allocated space to hold values for')
                        end
                        
                        for trajCtr = 1:size(continuousFrameRuns,2)
                            continuousframes = continuousFrameRuns{trajCtr};
                            wormtraj_xcoords = NaN(numel(continuousframes),size(worm_xcoords,2));
                            wormtraj_ycoords = NaN(numel(continuousframes),size(worm_ycoords,2));
                            for trajRunFrameCtr = 1:numel(continuousframes)
                                wormtraj_xcoords(trajRunFrameCtr,:) = worm_xcoords(wormCatLogInd...
                                    & trajData.frame_number == continuousframes(trajRunFrameCtr),:);
                                wormtraj_ycoords(trajRunFrameCtr,:) = worm_ycoords(wormCatLogInd...
                                    & trajData.frame_number == continuousframes(trajRunFrameCtr),:);
                            end
                            
                            % filter for maximum traj length
                            maxTrajFrameNum = round(maxTrajDuration*frameRate);
                            if size(wormtraj_xcoords,1)> maxTrajFrameNum
                                if strcmp(wormcats{wormcatCtr},'leaveCluster')
                                    % always start leave cluster traj from the moment the worm exits cluster
                                    firstFrame = 1;
                                elseif strcmp(wormcats{wormcatCtr},'loneWorm')
                                    % randomly sample the start of lone worm traj
                                    firstFrame = randi(size(wormtraj_xcoords,1)-maxTrajFrameNum,1);
                                end
                                % truncate the traj at maximum length
                                lastFrame = firstFrame + maxTrajFrameNum;
                                wormtraj_xcoords = wormtraj_xcoords(firstFrame:lastFrame,:);
                                wormtraj_ycoords = wormtraj_ycoords(firstFrame:lastFrame,:);
                            end
                            set(0,'CurrentFigure',fullTraj)
                            if strcmp(wormcats{wormcatCtr},'leaveCluster')
                                plot(wormtraj_xcoords,wormtraj_ycoords,'r')
                            elseif strcmp(wormcats{wormcatCtr},'loneWorm')
                                plot(wormtraj_xcoords,wormtraj_ycoords,'g')
                            end
                        end
                    end
                end
                set(0,'CurrentFigure',fullTraj)
                title([filename(end-31:end-18) ', worm' num2str(uniqueWorms(wormCtr)) '\_trajectory'])
                legend('full traj','leave cluster', 'lone worm')
                axis equal
                xlabel('microns')
                ylabel('microns')
                figurename = (['figures/turns/redManualTrajFull/redManualTrajFull_' phase '_' filename(end-31:end-18) '_worm' num2str(uniqueWorms(wormCtr))]);
                exportfig(fullTraj,[figurename '.eps'],exportOptions)
                system(['epstopdf ' figurename '.eps']);
                system(['rm ' figurename '.eps']);
            end
        end
    end
end