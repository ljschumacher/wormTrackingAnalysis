clear
close all

addpath('auxiliary/')

% Script is based on visualizeRedManualTraj.m, except it lets you specify
% first and last frames for plotting specific catergories of worms. This is
% to help visualize trajectories of interest and narrow down frame ranges
% (which is useful for observing the event from actual data). 

% Script uses manually joined bodywall muscle marker trajectory data and
% goes worm by worm to plot its full trajectory (for the specified phase)
% in black, and then categorised trajectories (leave cluster vs. lone worm)
% using specified coloring schemes. The options are to color according to
% movement direction, according to signed speed, or simply based on
% category (if neither of the other option is selected).

%% set parameters
phase = 'partial'; % 'partial' to indicate that for this particular script, frames of interest are specified by hand
marker = 'bodywall'; % 'bodywall'
strains = {'npr1'}; % {'npr1','N2'}
wormnums = {'40'};% {'40'};
wormcats = {'leaveCluster','loneWorm'}; %'leaveCluster','loneWorm'
preExitDuration = 2; % only applied if colorSpeed is true: duration (in seconds) before a worm exits a cluster to be included in the leave cluster analysis
postExitDuration = 5; % duration (in seconds) after a worm exits a cluster to be included in the leave cluster analysis
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
colorDirection = false;
colorSpeed = true;
saveResults = false;

currentFile = 5; %
firstFrame = 18800; %
lastFrame = 19000; %

maxBlobSize_r = 2.5e5;
minSkelLength_r = 850;
maxSkelLength_r = 1500;
minNeighbrDist = 2000;
inClusterNeighbourNum = 3;

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        % load file list
        [~,filenames,~] = xlsread(['datalists/' strains{strainCtr} '_' wormnum '_r_list.xlsx'],1,'A1:E15','basic');

        %% go through individual movies
        for fileCtr = currentFile
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
            worm_xcoords = mean(worm_xcoords(:,1:8),2);
            worm_ycoords = mean(worm_ycoords(:,1:8),2);
            
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
            % apply frame restriction            
            phaseFrameLogInd = trajData.frame_number <= lastFrame & trajData.frame_number >= firstFrame;
            trajData.filtered(~phaseFrameLogInd) = false;
            % find worms that have just left a cluster vs lone worms
            [leaveClusterLogInd, loneWormLogInd,~,~] = findWormCategory(filename,inClusterNeighbourNum,minNeighbrDist,postExitDuration);
            if colorSpeed
                % for colorSpeed, turn on a few seconds prior to cluster exit for leaveClusterLogInd
                [~,leaveClusterLogInd,~,~] = findClusterEntryExitX(filename,inClusterNeighbourNum,minNeighbrDist,preExitDuration,postExitDuration);
            end
            
            %% get movement direction
            if colorDirection | colorSpeed
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
                
                midbodySpeedSigned = midbodySpeedSigned';
                fwdLogInd = midbodySpeedSigned>0;
                revLogInd = midbodySpeedSigned<0;
                nanLogInd = isnan(midbodySpeedSigned);
                
            end
            
            %% plot worm full vs. categorised traj's
            
            % go worm by worm and plot full traj
            uniqueWorms = unique(features.worm_index);
            for wormCtr = 1:numel(uniqueWorms)
                uniqueWorm = uniqueWorms(wormCtr);
                % use unique worm
                wormLogInd = trajData.worm_index_manual==uniqueWorms(wormCtr);
                % apply phase restriction
                phaseFrameLogInd = trajData.frame_number <= lastFrame & trajData.frame_number >= firstFrame;
                wormLogInd(~phaseFrameLogInd) = false;
                % get xy coordinates
                wormtraj_xcoords = worm_xcoords(wormLogInd);
                wormtraj_ycoords = worm_ycoords(wormLogInd);
                % plot full trajectory in black
                fullTraj = figure; hold on
                plot(wormtraj_xcoords,wormtraj_ycoords,'k')
                % now plot individually categorised traj on top of the full traj
                wormLogInd = trajData.worm_index_manual==uniqueWorms(wormCtr) & trajData.filtered;
                % go through each worm category
                for wormcatCtr = 1:length(wormcats)
                    wormCatLogInd = wormLogInd & eval([wormcats{wormcatCtr} 'LogInd']);
                    % plot xy coordinates for categorised worms
                    wormtraj_xcoords = worm_xcoords(wormCatLogInd);
                    wormtraj_ycoords = worm_ycoords(wormCatLogInd);
                    if ~(colorDirection | colorSpeed)
                        set(0,'CurrentFigure',fullTraj)
                        if strcmp(wormcats{wormcatCtr},'leaveCluster')
                            plot(wormtraj_xcoords,wormtraj_ycoords,'r.','MarkerSize',10)
                        elseif strcmp(wormcats{wormcatCtr},'loneWorm')
                            plot(wormtraj_xcoords,wormtraj_ycoords,'b.','MarkerSize',10)
                        end
                    elseif colorSpeed
                        colors = midbodySpeedSigned(wormCatLogInd);
                        set(0,'CurrentFigure',fullTraj)
                        % plot leaveCluster worm only
                        if strcmp(wormcats{wormcatCtr},'leaveCluster')
                            scatter(wormtraj_xcoords,wormtraj_ycoords,8,colors,'filled')
                        end
                        load('divergentColorMap.txt')
                        colormap(divergentColorMap)
                        caxis([-500 500])
                        colorbar
                    elseif colorDirection
                        % optional: color trajectory by direction of movement
                        wormtraj_xcoords_fwd = worm_xcoords(wormCatLogInd & fwdLogInd);
                        wormtraj_ycoords_fwd = worm_ycoords(wormCatLogInd & fwdLogInd);
                        wormtraj_xcoords_rev = worm_xcoords(wormCatLogInd & revLogInd);
                        wormtraj_ycoords_rev = worm_ycoords(wormCatLogInd & revLogInd);
                        wormtraj_xcoords_nan = worm_xcoords(wormCatLogInd & nanLogInd);
                        wormtraj_ycoords_nan = worm_ycoords(wormCatLogInd & nanLogInd);
                        set(0,'CurrentFigure',fullTraj)
                        if strcmp(wormcats{wormcatCtr},'leaveCluster')
                            plot(wormtraj_xcoords_fwd,wormtraj_ycoords_fwd,'r.','MarkerSize',10)
                            plot(wormtraj_xcoords_rev,wormtraj_ycoords_rev,'c.','MarkerSize',10)
                            plot(wormtraj_xcoords_nan,wormtraj_ycoords_nan,'y.','MarkerSize',10)
                        elseif strcmp(wormcats{wormcatCtr},'loneWorm')
                            plot(wormtraj_xcoords_fwd,wormtraj_ycoords_fwd,'m.','MarkerSize',10)
                            plot(wormtraj_xcoords_rev,wormtraj_ycoords_rev,'b.','MarkerSize',10)
                            plot(wormtraj_xcoords_nan,wormtraj_ycoords_nan,'y.','MarkerSize',10)
                        end 
                    end
                end
                % format and save figures
                set(0,'CurrentFigure',fullTraj)
                title([filename(end-31:end-18) ', worm' num2str(uniqueWorms(wormCtr)) '\_trajectory'])
                axis equal
                xlabel('microns')
                ylabel('microns')
                xlim([0 12000])
                ylim([0 12000])
                if colorDirection
                    figurename = (['figures/redManualTrajFull/redManualTrajFullDir_' phase '_' filename(end-31:end-18) '_worm' num2str(uniqueWorms(wormCtr))]);
                elseif colorSpeed
                    figurename = (['figures/redManualTrajFull/redManualTrajFullSpeed_' phase '_' filename(end-31:end-18) '_worm' num2str(uniqueWorms(wormCtr))]);
                else
                    figurename = (['figures/redManualTrajFull/redManualTrajFull_' phase '_' filename(end-31:end-18) '_worm' num2str(uniqueWorms(wormCtr))]);
                end
                if saveResults
                exportfig(fullTraj,[figurename '.eps'],exportOptions)
                system(['epstopdf ' figurename '.eps']);
                system(['rm ' figurename '.eps']);
                end
            end
        end
    end
end