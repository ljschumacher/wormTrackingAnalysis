% plot features of red worms (from the second dataset) that may be useful
% for identifying omega turns

clear
close all

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

%% set parameters
phase = 'fullMovie'; % 'fullMovie', 'joining', or 'sweeping'.
strains = {'npr1','N2'}; % {'npr1','N2'}
wormnums = {'40'};% {'40','HD'};

intensityThresholds = containers.Map({'40','HD','1W'},{60, 40, 100});
maxBlobSize = 2.5e5;
minSkelLength = 850;
maxSkelLength = 1500;
minNeighbrDist = 2000;
inClusterNeighbourNum = 3;
pixelsize = 100/19.5; % 100 microns are 19.5 pixels

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        % load data
        [phaseFrames,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_r_list.xlsx'],1,'A1:E15','basic');
        phaseFrames = phaseFrames-1; % to correct for python indexing at 0
        numFiles = length(filenames);
        % create cell arrays to hold individual movie values to be pooled
        maxPathNum = 200; % assume maximum number of paths per file is as set. Script checks for this later and gives warning if this is not enough
%         intHeadAngles_all = cell(1,maxPathNum,numFiles);
%         intHeadAngles_leaveCluster = cell(1,maxPathNum,numFiles);
%         intHeadAngles_loneWorms = cell(1,maxPathNum,numFiles);
        %         headAngularSpeed_leaveCluster = cell(numFiles,1);
        %         headAngularSpeed_loneWorms = cell(numFiles,1);
        %         omegaTurns_leaveCluster = cell(numFiles,1);
        %         omegaTurns_loneWorms = cell(numFiles,1);
        intHeadAngleByChunk_allFiles = cell(1,numFiles);
        intHeadAnglesFig = figure;
        %headAngularSpeedFig = figure;
        %omegaTurnsFig = figure;
        leaveClusterWormCount = 0;
        loneWormCount = 0;
        %% go through individual movies
        for fileCtr = 1:numFiles
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton');
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            features = h5read(strrep(filename,'skeletons','features'),'/features_timeseries');
            if strcmp(phase, 'fullMovie')
                firstFrame = 0;
                lastFrame = phaseFrames(fileCtr,4);
            elseif strcmp(phase,'joining')
                firstFrame = phaseFrames(fileCtr,1);
                lastFrame = phaseFrames(fileCtr,2);
            elseif strcmp(phase,'sweeping')
                firstFrame = phaseFrames(fileCtr,3);
                lastFrame = phaseFrames(fileCtr,4);
            end
            %% filter worms by various criteria
            % filter red by blob size and intensity
            if contains(filename,'55')||contains(filename,'54')
                intensityThreshold = 80;
            else
                intensityThreshold = 40;
            end
            trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                intensityThreshold,maxBlobSize);
            % filter red by skeleton length
            trajData.filtered = trajData.filtered&logical(trajData.is_good_skel)...
                &filterSkelLength(skelData,pixelsize,minSkelLength,maxSkelLength);
            % apply phase restriction
            phaseFrameLogInd = trajData.frame_number <= lastFrame & trajData.frame_number >= firstFrame;
            trajData.filtered(~phaseFrameLogInd) = false;
            features.filtered = ismember(features.skeleton_id+1,find(trajData.filtered)); % use trajData.filtered to filter out unwa
            % find worms that have just left a cluster
            min_neighbr_dist = h5read(filename,'/min_neighbr_dist');
            num_close_neighbrs = h5read(filename,'/num_close_neighbrs');
            neighbr_dist = h5read(filename,'/neighbr_distances');
            inClusterLogInd = num_close_neighbrs>=inClusterNeighbourNum;
            leaveClusterLogInd = vertcat(false,inClusterLogInd(1:end-1)&~inClusterLogInd(2:end)); % find worm-frames where inCluster changes from true to false
            leaveClusterFrameStart = find(leaveClusterLogInd);
            leaveClusterFrameEnd = leaveClusterFrameStart+5*frameRate; %retain data for 5 seconds after a worm exits cluster
            leaveClusterFrameEnd = leaveClusterFrameEnd(leaveClusterFrameEnd<=numel(leaveClusterLogInd)); % exclude movie segments with frames beyond highest frame number
            leaveClusterFrameStart = leaveClusterFrameStart(1:numel(leaveClusterFrameEnd));
            for exitCtr = 1:numel(leaveClusterFrameStart)
                leaveClusterLogInd(leaveClusterFrameStart(exitCtr):leaveClusterFrameEnd(exitCtr))=true;
            end
            %the following line should do the same as the loop above, if dimensions of
            %terms added are compatible (eg row + column vector)
            %leaveClusterLogInd(unique(leaveClusterFrameStart + 0:5*frameRate)) = true;
            leaveClusterLogInd(inClusterLogInd)=false; % exclude when worms move back into a cluster
            loneWormLogInd = min_neighbr_dist>=minNeighbrDist;
            leaveClusterLogInd(loneWormLogInd)=false; % exclude worms that have become lone worm
            leaveClusterLogInd = ismember(features.skeleton_id+1,find(trajData.filtered & leaveClusterLogInd)); % make clusterClusterLogInd the same size as features.filtered
            % find lone worms
            loneWormLogInd = ismember(features.skeleton_id+1,find(trajData.filtered & loneWormLogInd));
            %% calculate or extract desired feature values
            % calculte integrated angles
            uniquePaths = unique(features.worm_index);
            intHeadAngleByChunk_allPaths = cell(1,numel(uniquePaths));
            for pathCtr = 1:numel(uniquePaths)
                path = uniquePaths(pathCtr);
                pathLogInd = ismember(features.worm_index, path);
                headAnglesAllLogInd = pathLogInd & features.filtered;
%                 headAngles_all = features.head_tip_motion_direction(headAnglesAllLogInd);
%                 intHeadAngles_all{1,pathCtr,fileCtr} = abs(cumsum(headAngles_all));
                headAnglesLeaveClusterLogInd = headAnglesAllLogInd & leaveClusterLogInd;
                headAngles_leaveCluster = features.head_tip_motion_direction(headAnglesLeaveClusterLogInd);
                if size(headAngles_leaveCluster,1)> 10*frameRate % only concerned about turn events within 10 seconds of exiting cluster
                    headAngles_leaveCluster = headAngles_leaveCluster(1:10*frameRate);
                end
                headAngles_leaveCluster(isnan(headAngles_leaveCluster))=[];
                if ~isempty(headAngles_leaveCluster)
                    headAngles_leaveClusterChunks = vertcat(1, find(headAngles_leaveCluster(1:end-1).*...
                        headAngles_leaveCluster(2:end)<0), numel(headAngles_leaveCluster));
                    % generate a list of positions of headAngles_leaveCluster vector where signs change
                    numChunk = length(headAngles_leaveClusterChunks)-1;
                    intHeadAngleByChunk_singlePath = NaN(1,numChunk); % create variable to hold integrated values from each continuous turn pattern
                    for chunkCtr = 1:numChunk
                        chunkIntHeadAngle = abs(sum(headAngles_leaveCluster(headAngles_leaveClusterChunks(chunkCtr)):...
                            headAngles_leaveCluster(headAngles_leaveClusterChunks(chunkCtr+1))));
                        intHeadAngleByChunk_singlePath(chunkCtr) = chunkIntHeadAngle;
                    end
                    intHeadAngleByChunk_allPaths{pathCtr}=intHeadAngleByChunk_singlePath;
                end
%                 headAnglesLoneWormLogInd = headAnglesAllLogInd & loneWormLogInd;
%                 headAngles_loneWorm = features.head_tip_motion_direction(headAnglesLoneWormLogInd);
%                 intHeadAngles_loneWorms{1,pathCtr,fileCtr}= abs(cumsum(headAngles_loneWorm));
%                 if ~isempty(intHeadAngles_all{:,maxPathNum,fileCtr})
%                     warning(['some worm paths may be missed off from integrated head angle calculations in' filename])
%                 end
%             end
            intHeadAngleByChunk_allPathsCat = horzcat(intHeadAngleByChunk_allPaths{:});
            intHeadAngleByChunk_allFiles{fileCtr} = intHeadAngleByChunk_allPathsCat;
            % write additional turn features into cell array for pooling across movies later
            %headAngularSpeed_leaveCluster{fileCtr} = abs(features.head_motion_direction(features.filtered&leaveClusterLogInd));
            %$$$$$$$$$$$$ omegaTurns_leaveCluster{fileCtr} = abs(features.path_curvature(features.filtered&leaveClusterLogInd));
            %headAngularSpeed_loneWorms{fileCtr} = abs(features.head_motion_direction(features.filtered&loneWormLogInd));
            %$$$$$$$$$$$$ omegaTurns_loneWorms{fileCtr} = abs(features.path_curvature(features.filtered&loneWormLogInd));
            % add number of worms for total n number across movies
            %leaveClusterWormCount = leaveClusterWormCount + numel(unique(features.worm_index(leaveClusterLogInd)));
            %loneWormCount = loneWormCount + numel(unique(features.worm_index(loneWormLogInd)));
            end
        end
        intHeadAngleByChunk_allFiles = horzcat(intHeadAngleByChunk_allFiles{:})';
        intHeadAngleByChunk_allFiles = intHeadAngleByChunk_allFiles(intHeadAngleByChunk_allFiles>0);
    end
end
        %% pool data from all files belonging to the same strain and worm density
%         %
%         intHeadAngles_all = intHeadAngles_all(~cellfun('isempty',intHeadAngles_all));
%         [longestPath_all,~] = max(cellfun('size',intHeadAngles_all,1));
%         intHeadAngles_all_C2M = NaN(size(intHeadAngles_all,1),longestPath_all); % turn cell array into matrix
%         for pathPlotCtr = 1:size(intHeadAngles_all,1)
%             intHeadAngles_all_C2M(pathPlotCtr,1:length(intHeadAngles_all{pathPlotCtr}')) = intHeadAngles_all{pathPlotCtr}';
%         end
%         %
%         intHeadAngles_leaveCluster = intHeadAngles_leaveCluster...
%             (~cellfun('isempty',intHeadAngles_leaveCluster));
%         [longestPath_leaveCluster,~] = max(cellfun('size',intHeadAngles_leaveCluster,1));
%         intHeadAngles_leaveCluster_C2M = NaN(size(intHeadAngles_leaveCluster,1),...
%             longestPath_leaveCluster); % turn cell array into matrix
%         for pathPlotCtr = 1:size(intHeadAngles_leaveCluster,1)
%             intHeadAngles_leaveCluster_C2M(pathPlotCtr,1:length(intHeadAngles_leaveCluster{pathPlotCtr}'))...
%                 = intHeadAngles_leaveCluster{pathPlotCtr}';
%         end
%         %
%         intHeadAngles_loneWorms = intHeadAngles_loneWorms...
%             (~cellfun('isempty',intHeadAngles_loneWorms));
%         [longestPath_loneWorms,~] = max(cellfun('size',intHeadAngles_loneWorms,1));
%         intHeadAngles_loneWorms_C2M = NaN(size(intHeadAngles_loneWorms,1),...
%             longestPath_loneWorms); % turn cell array into matrix
%         for pathPlotCtr = 1:size(intHeadAngles_loneWorms,1)
%             intHeadAngles_loneWorms_C2M(pathPlotCtr,1:length(intHeadAngles_loneWorms{pathPlotCtr}'))...
%                 = intHeadAngles_loneWorms{pathPlotCtr}';
%         end
        %
        %         headAngularSpeed_leaveCluster = vertcat(headAngularSpeed_leaveCluster{:});
        %         headAngularSpeed_loneWorms = vertcat(headAngularSpeed_loneWorms{:});
        %omegaTurns_leaveCluster = vertcat(omegaTurns_leaveCluster{:});
        %omegaTurns_loneWorms = vertcat(omegaTurns_loneWorms{:});
        %% plot data
%         % intHeadAngles figure
%         set(0,'CurrentFigure',intHeadAnglesFig); hold on
%         shadedErrorBar(1:longestPath_all,nanmean(intHeadAngles_all_C2M,1),...
%             nanstd(intHeadAngles_all_C2M,1)./sqrt(sum(~isnan(intHeadAngles_all_C2M))),'r',1)
%         shadedErrorBar(1:longestPath_leaveCluster,nanmean(intHeadAngles_leaveCluster_C2M,1),...
%             nanstd(intHeadAngles_leaveCluster_C2M,1)./sqrt(sum(~isnan(intHeadAngles_leaveCluster_C2M))),'b',1)
%         shadedErrorBar(1:longestPath_loneWorms,nanmean(intHeadAngles_loneWorms_C2M,1),...
%             nanstd(intHeadAngles_loneWorms_C2M,1)./sqrt(sum(~isnan(intHeadAngles_loneWorms_C2M))),'k',1)
%         %legend('all','','','','leave cluster','','','','lone worms','','','')
%         xlabel('number of frames')
%         %xlim([0 5400])
%         ylabel('mean integrated head angles')
%         set(intHeadAnglesFig,'PaperUnits','centimeters')
%         figurename = ['figures/turns/intHeadAngles_' strains{strainCtr} '_' wormnums{numCtr} '_' phase];
%         exportfig(intHeadAnglesFig,[figurename '.eps'],exportOptions)
%         system(['epstopdf ' figurename '.eps']);
%         system(['rm ' figurename '.eps']);
        
        % headAngularSpeed figure
        %         set(0,'CurrentFigure',headAngularSpeedFig)
        %         histogram(headAngularSpeed_leaveCluster,'Normalization','pdf','DisplayStyle','stairs')
        %         hold on
        %         histogram(headAngularSpeed_loneWorms,'Normalization','pdf','DisplayStyle','stairs')
        %         leaveClusterLegend = ['leave cluster, n = ' num2str(leaveClusterWormCount)];
        %         loneWormLegend = ['lone worm, n = ' num2str(loneWormCount)];
        %         legend(leaveClusterLegend, loneWormLegend)
        %         title([strains{strainCtr} '\_' wormnums{numCtr} '\_headAngularSpeed'],'FontWeight','normal')
        %         xlabel('head angular speed (degree/s)')
        %         ylabel('probability')
        %         xlim([0 8])
        %         ylim([0 2])
        %         set(headAngularSpeedFig,'PaperUnits','centimeters')
        %         figurename = ['figures/turns/headAngularSpeed_' strains{strainCtr} '_' wormnums{numCtr} '_' phase '_CL'];
        %         exportfig(headAngularSpeedFig,[figurename '.eps'],exportOptions)
        %         system(['epstopdf ' figurename '.eps']);
        %         system(['rm ' figurename '.eps']);
        
        % omegaTurns figure
        %         set(0,'CurrentFigure',omegaTurnsFig)
        %         histogram(omegaTurns_leaveCluster,'Normalization','pdf','DisplayStyle','stairs')
        %         hold on
        %         histogram(omegaTurns_loneWorms,'Normalization','pdf','DisplayStyle','stairs')
        %         legend(leaveClusterLegend, loneWormLegend)
        %         title([strains{strainCtr} '\_' wormnums{numCtr} '\_omegaTurns'],'FontWeight','normal')
        %         xlabel('path curvature (radians/microns)')
        %         ylabel('probability')
        %         xlim([0 0.5])
        %         ylim([0 60])
        %         set(omegaTurnsFig,'PaperUnits','centimeters')
        %         figurename = ['figures/turns/omegaTurns_' strains{strainCtr} '_' wormnums{numCtr} '_' phase '_CL'];
        %         exportfig(omegaTurnsFig,[figurename '.eps'],exportOptions)
        %         system(['epstopdf ' figurename '.eps']);
        %         system(['rm ' figurename '.eps']);
%     end
% end