clear
close all

% Script plots midbody speed for worms entering or leaving a cluster using the manually joined red trajectories, with
% the point of entry/exit aligned.

%% set parameters
phase = 'joining'; % 'fullMovie', 'joining', or 'sweeping'.
strains = {'npr1'}; % {'npr1','N2'}
wormnums = {'40'};% {'40'};
preExitDuration = 10; % only applied if colorSpeed is true: duration (in seconds) before a worm exits a cluster to be included in the leave cluster analysis
postExitDuration = 10; % duration (in seconds) after a worm exits a cluster to be included in the leave cluster analysis
minInOutClusterFrameNum = 5;
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
saveResults = true;
enforceInClusterAfterEntryBeforeExit = false;
smoothWindow = 5;

maxBlobSize_r = 2.5e5;
minSkelLength_r = 850;
maxSkelLength_r = 1500;
minNeighbrDist = 2000;
inClusterNeighbourNum = 3;

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',30,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',20,...
    'LineWidth',1);

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        % load file list
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_r_list.xlsx'],1,'A1:E15','basic');
        numFiles = length(filenames);
        
        %% go through individual movies
        for fileCtr = 1:numFiles
            %% load data
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton'); % in pixels
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            features = h5read(strrep(filename,'skeletons','feat_manual'),'/features_timeseries');
            
            %% generate time series (x), setting time to be zero at the point of entry
            timeSeries = [-preExitDuration*frameRate:postExitDuration*frameRate];
            
            %% calculate midbody signed speed (from reversalAnalysisBodyWall.m)
            midbodyIndcs = 19:33;
            % centroids of midbody skeleton
            midbody_x = mean(squeeze(skelData(1,midbodyIndcs,:)))*pixelsize;
            midbody_y = mean(squeeze(skelData(2,midbodyIndcs,:)))*pixelsize;
            % change in centroid position over time (issue of worm shifts?)
            dmidbody_xdt = gradient(midbody_x)*frameRate;
            dmidbody_ydt = gradient(midbody_y)*frameRate;
            % midbody speed and velocity
            dFramedt = gradient(double(trajData.frame_number))';
            midbodySpeed = sqrt(dmidbody_xdt.^2 + dmidbody_ydt.^2)./dFramedt;
            midbodyVelocity = [dmidbody_xdt; dmidbody_ydt]./dFramedt;
            % direction of segments pointing along midbody
            [~, dmidbody_yds] = gradient(squeeze(skelData(2,midbodyIndcs,:)),-1);
            [~, dmidbody_xds] = gradient(squeeze(skelData(1,midbodyIndcs,:)),-1);
            % sign speed based on relative orientation of velocity to midbody
            midbodySpeedSigned = getSignedSpeed(midbodyVelocity,[mean(dmidbody_xds); mean(dmidbody_yds)]);
            % ignore first and last frames of each worm's track
            wormChangeIndcs = gradient(double(trajData.worm_index_manual))~=0;
            midbodySpeedSigned(wormChangeIndcs)=NaN;        
            
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
            % get logical indices for worms that have entered or left a cluster
            min_neighbr_dist = h5read(filename,'/min_neighbr_dist');
            num_close_neighbrs = h5read(filename,'/num_close_neighbrs');
            trajData = h5read(filename,'/trajectories_data');
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            inClusterLogInd = num_close_neighbrs>=inClusterNeighbourNum;
            enterClusterStartLogInd = vertcat(false,~inClusterLogInd(1:end-1)&inClusterLogInd(2:end));
            exitClusterStartLogInd = vertcat(false,inClusterLogInd(1:end-1)&~inClusterLogInd(2:end));
            
            %% plot entry time course for midbody signed speed      
            % get indices for cluster entry points
            enterClusterPoint = find(enterClusterStartLogInd);  
            entrySpeeds = NaN(numel(enterClusterPoint),preExitDuration*frameRate +1+ postExitDuration*frameRate);
            % loop through each entry
            for entryCtr = 1:numel(enterClusterPoint)
                thisEntryIdx = enterClusterPoint(entryCtr);
                % get the worm index
                wormIndex = trajData.worm_index_manual(thisEntryIdx);
                % check that the same worm stays in cluster for minumum number of frames
                if thisEntryIdx+minInOutClusterFrameNum <= length(trajData.frame_number) &&...
                        nnz(inClusterLogInd(thisEntryIdx:(thisEntryIdx+minInOutClusterFrameNum))) == 1+minInOutClusterFrameNum...
                        && nnz(trajData.worm_index_manual(thisEntryIdx:(thisEntryIdx+minInOutClusterFrameNum)) == wormIndex) == 1+minInOutClusterFrameNum
                    % expand for a specified number of entries before and after the entry point
                    thisEntryStartIdx = thisEntryIdx-preExitDuration*frameRate;
                    if thisEntryStartIdx <= 0
                        beforeStartFrameNum = 1-thisEntryStartIdx; % take note of omitted frames for alignment purposes
                        thisEntryStartIdx = 1; % exclude entries below index 0
                    end
                    thisEntryEndIdx = thisEntryIdx + postExitDuration*frameRate;
                    if thisEntryEndIdx > length(trajData.frame_number)
                        afterEndFrameNum = thisEntryEndIdx - length(trajData.frame_number); % take note of omitted frames for alignment purposes
                        thisEntryEndIdx = length(trajData.frame_number); % exclude entries above highest index
                    end
                    % generate logical index for this entry
                    thisEntryLogInd = false(1,length(midbodySpeedSigned));
                    thisEntryLogInd(thisEntryStartIdx:thisEntryEndIdx) = true;
                    % get midbody signed speeds (y) for this entry
                    thisEntrySpeeds = NaN(1,size(entrySpeeds,2));
                    if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                        thisEntrySpeeds(beforeStartFrameNum+1:end) = midbodySpeedSigned(thisEntryLogInd);
                    elseif exist('afterEndFrameNum','var')
                        thisEntrySpeeds(1:(end-afterEndFrameNum)) = midbodySpeedSigned(thisEntryLogInd);
                    else
                        thisEntrySpeeds = midbodySpeedSigned(thisEntryLogInd);
                    end
                    % exclude entries representing a different worm
                    thisEntrySpeeds(trajData.worm_index_manual(thisEntryLogInd) ~= wormIndex) = NaN;
                    % optional: excludes entries representing non-inCluster worms post-entry (inclusive of entry point)
                    if enforceInClusterAfterEntryBeforeExit
                        % create speed vector the same length as the full index
                        thisEntrySpeedsFullLength = NaN(1,length(midbodySpeedSigned));
                        if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                             thisEntrySpeedsFullLength(thisEntryStartIdx:thisEntryEndIdx) = thisEntrySpeeds(beforeStartFrameNum+1:end);
                        elseif exist('afterEndFrameNum','var') 
                            thisEntrySpeedsFullLength(thisEntryStartIdx:thisEntryEndIdx) = thisEntrySpeeds(1:(end-afterEndFrameNum));
                        else
                            thisEntrySpeedsFullLength(thisEntryStartIdx:thisEntryEndIdx) = thisEntrySpeeds;
                        end
                        % generate logical index for post-entry indices
                        thisEntryPostEntryLogInd = false(1,length(midbodySpeedSigned));
                        thisEntryPostEntryLogInd(thisEntryIdx:thisEntryEndIdx) = true;
                        % eliminate speeds for non-inCluster worms post-entry
                        thisEntrySpeedsFullLength(thisEntryPostEntryLogInd & ~ inClusterLogInd') = NaN;
                        % shorten full length speeds back to the window of interest
                        if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                            thisEntrySpeeds(beforeStartFrameNum+1:end) = thisEntrySpeedsFullLength(thisEntryStartIdx:thisEntryEndIdx);
                        elseif exist('afterEndFrameNum','var') 
                            thisEntrySpeeds(1:(end-afterEndFrameNum)) = thisEntrySpeedsFullLength(thisEntryStartIdx:thisEntryEndIdx);
                        else
                            thisEntrySpeeds = thisEntrySpeedsFullLength(thisEntryStartIdx:thisEntryEndIdx);
                        end
                    end
                    % add data to speed matrix
                    entrySpeeds(entryCtr,:) = thisEntrySpeeds;
                end
                % clear variables
                clear beforeStartFrameNum
                clear afterEndFrameNum
            end
            % set maximum speed
            entrySpeeds(entrySpeeds>1500) = NaN;
            % smooth speeds
            smoothEntrySpeeds = smoothdata(entrySpeeds,2,'movmean',smoothWindow,'includenan');
            % set maximum speed
            smoothEntrySpeeds(smoothEntrySpeeds>1500) = NaN;

            %% plot exit time course for midbody signed speed
            % get indices for cluster exit points
            exitClusterPoint = find(exitClusterStartLogInd);
            exitSpeeds = NaN(numel(exitClusterPoint),preExitDuration*frameRate +1+ postExitDuration*frameRate);
            % loop through each exit
            for exitCtr = 1:numel(exitClusterPoint)
                thisExitIdx = exitClusterPoint(exitCtr);
                % get the worm index
                wormIndex = trajData.worm_index_manual(thisExitIdx);
                % check that the same worm stays out of cluster for minumum number of frames
                if thisExitIdx+minInOutClusterFrameNum <= length(trajData.frame_number) &&...
                    nnz(inClusterLogInd(thisExitIdx:(thisExitIdx+minInOutClusterFrameNum))) == 0 ...
                        && nnz(trajData.worm_index_manual(thisExitIdx:(thisExitIdx+minInOutClusterFrameNum)) == wormIndex) == 1+minInOutClusterFrameNum
                    % expand for a specified number of entries before and after the exit point
                    thisExitStartIdx = thisExitIdx-preExitDuration*frameRate;
                    if thisExitStartIdx <= 0
                        beforeStartFrameNum = 1-thisExitStartIdx; % take note of omitted frames for alignment purposes
                        thisExitStartIdx = 1; % exclude entries below index 0
                    end
                    thisExitEndIdx = thisExitIdx + postExitDuration*frameRate;
                    if thisExitEndIdx > length(trajData.frame_number)
                        afterEndFrameNum = thisExitEndIdx - length(trajData.frame_number); % take note of omitted frames for alignment purposes
                        thisExitEndIdx = length(trajData.frame_number); % exclude entries above highest index
                    end
                    % generate logical index for this exit
                    thisExitLogInd = false(1,length(midbodySpeedSigned));
                    thisExitLogInd(thisExitStartIdx:thisExitEndIdx) = true;
                    % get midbody signed speeds (y) for this exit
                    thisExitSpeeds = NaN(1,size(exitSpeeds,2));
                    if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                        thisExitSpeeds(beforeStartFrameNum+1:end) = midbodySpeedSigned(thisExitLogInd);
                    elseif exist('afterEndFrameNum','var')
                        thisExitSpeeds(1:(end-afterEndFrameNum)) = midbodySpeedSigned(thisExitLogInd);
                    else
                        thisExitSpeeds = midbodySpeedSigned(thisExitLogInd);
                    end
                    % exclude entries representing a different worm
                    thisExitSpeeds(trajData.worm_index_manual(thisExitLogInd) ~= wormIndex) = NaN;
                    % optional: excludes entries representing non-inCluster worms post-entry (inclusive of entry point)
                    if enforceInClusterAfterEntryBeforeExit
                        % create speed vector the same length as the full index
                        thisExitSpeedsFullLength = NaN(1,length(midbodySpeedSigned));
                        if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                            thisExitSpeedsFullLength(thisExitStartIdx:thisExitEndIdx) = thisExitSpeeds(beforeStartFrameNum+1:end);
                        elseif exist('afterEndFrameNum','var')
                            thisExitSpeedsFullLength(thisExitStartIdx:thisExitEndIdx) = thisExitSpeeds(1:(end-afterEndFrameNum));
                        else
                            thisExitSpeedsFullLength(thisExitStartIdx:thisExitEndIdx) = thisExitSpeeds;
                        end
                        % generate logical index for post-entry indices
                        thisExitPreExitLogInd = false(1,length(midbodySpeedSigned));
                        thisExitPreExitLogInd(thisExitStartIdx:(thisExitIdx-1)) = true;
                        % eliminate speeds for non-inCluster worms post-entry
                        thisExitSpeedsFullLength(thisExitPreExitLogInd & ~ inClusterLogInd') = NaN;
                        % shorten full length speeds back to the window of interest
                        if exist('beforeStartFrameNum','var') % this keeps the alignment of the entries in case they go below or above min/max index
                            thisExitSpeeds(beforeStartFrameNum+1:end) = thisExitSpeedsFullLength(thisExitStartIdx:thisExitEndIdx);
                        elseif exist('afterEndFrameNum','var')
                            thisExitSpeeds(1:(end-afterEndFrameNum)) = thisExitSpeedsFullLength(thisExitStartIdx:thisExitEndIdx);
                        else
                            thisExitSpeeds = thisExitSpeedsFullLength(thisExitStartIdx:thisExitEndIdx);
                        end
                    end
                    % add data to speed matrix
                    exitSpeeds(exitCtr,:) = thisExitSpeeds;
                end
                % clear variables
                clear beforeStartFrameNum
                clear afterEndFrameNum
            end
            % set maximum speed
            exitSpeeds(exitSpeeds>1500) = NaN;
            % smooth speeds
            smoothExitSpeeds = smoothdata(exitSpeeds,2,'movmean',smoothWindow,'includenan');
            % set maximum speed
            smoothExitSpeeds(smoothExitSpeeds>1500) = NaN;
            
            %% save data for pooling
            allEntrySpeeds{fileCtr} = entrySpeeds;
            allExitSpeeds{fileCtr} = exitSpeeds;
            allSmoothEntrySpeeds{fileCtr} = smoothEntrySpeeds;
            allSmoothExitSpeeds{fileCtr} = smoothExitSpeeds;

        end
        
        % pool data
        allEntrySpeeds = vertcat(allEntrySpeeds{:});
        allExitSpeeds = vertcat(allExitSpeeds{:});
        allSmoothEntrySpeeds = vertcat(allSmoothEntrySpeeds{:});
        allSmoothExitSpeeds = vertcat(allSmoothExitSpeeds{:});
        
        %% plotting and saving data
        % mean entry speed plot
        meanEntryExitSpeedsFig = figure; hold on
        set(0,'CurrentFigure',meanEntryExitSpeedsFig)
        plot(timeSeries,nanmean(allSmoothEntrySpeeds,1))
        plot(timeSeries,nanmean(allSmoothExitSpeeds,1))
        title(['mean cluster entry and exit speeds'])
        xlabel('frames')
        ylabel('speed(microns/s)')
        legend('entry','exit')
        if enforceInClusterAfterEntryBeforeExit
            figurename = (['figures/entryExitSpeeds/entryExitSpeedsMeanSmoothed_halfFiltered2_' phase]);
        else
            figurename = (['figures/entryExitSpeeds/entryExitSpeedsMeanSmoothed_' phase]);
        end
        if saveResults
            exportfig(meanEntryExitSpeedsFig,[figurename '.eps'],exportOptions)
            system(['epstopdf ' figurename '.eps']);
            system(['rm ' figurename '.eps']);
        end
        
        % errorbar entry plot
        figure;
        shadedErrorBar(timeSeries,nanmean(allSmoothEntrySpeeds,1),nanstd(allSmoothEntrySpeeds,1),'b');
        title(['mean cluster entry speeds'])
        xlabel('frames')
        ylabel('speed(microns/s)')
        ylim([-150 250])
        if enforceInClusterAfterEntryBeforeExit
            figurename = (['figures/entryExitSpeeds/entrySpeedsMeanErrorSmoothed_halfFiltered2_' phase]);
        else
            figurename = (['figures/entryExitSpeeds/entrySpeedsMeanErrorSmoothed_' phase]);
        end
        if saveResults
            entryErrorBarFig = gcf;
            exportfig(entryErrorBarFig,[figurename '.eps'],exportOptions)
            system(['epstopdf ' figurename '.eps']);
            system(['rm ' figurename '.eps']);
        end
        
        % errorbar exit plot
        figure;
        shadedErrorBar(timeSeries,nanmean(allSmoothExitSpeeds,1),nanstd(allSmoothExitSpeeds,1),'r');
        title(['mean cluster exit speeds'])
        xlabel('frames')
        ylabel('speed(microns/s)')
        ylim([-150 250])
        if enforceInClusterAfterEntryBeforeExit
            figurename = (['figures/entryExitSpeeds/exitSpeedsMeanErrorSmoothed_halfFiltered2_' phase]);
        else
            figurename = (['figures/entryExitSpeeds/exitSpeedsMeanErrorSmoothed_' phase]);
        end
        if saveResults
            exitErrorBarFig = gcf;
            exportfig(exitErrorBarFig,[figurename '.eps'],exportOptions)
            system(['epstopdf ' figurename '.eps']);
            system(['rm ' figurename '.eps']);
        end
        
        % save data
        if saveResults
            if enforceInClusterAfterEntryBeforeExit
                filename = './figures/entryExitSpeeds/allSmoothEntryExitSpeeds_halfFiltered2.mat';
            else
                filename = './figures/entryExitSpeeds/allSmoothEntryExitSpeeds.mat';
            end
            save(filename,'allSmoothEntrySpeeds','allSmoothExitSpeeds','allEntrySpeeds','allExitSpeeds')
        end
    end
end