
clear
close all

%% set parameters
phase = 'joining'; % 'fullMovie', 'joining', or 'sweeping'.
dataset = 2; % 1 or 2
marker = 'bodywall'; % 'pharynx' or 'bodywall'
strains = {'npr1'}; % {'npr1','N2'}
wormnums = {'40'};% {'40'};
wormcats = {'leaveCluster','loneWorm'}; %'leaveCluster','loneWorm'
postExitDuration = 5; % duration (in seconds) after a worm exits a cluster to be included in the leave cluster analysis
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
colorDirection = false;
colorSpeed = true;

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
        % create empty cell arrays to hold individual file values, so they can be pooled for a given strain/density combination
        for wormcatCtr = 1:length(wormcats)
            headAngTotal.(wormcats{wormcatCtr}) = cell(numFiles,1);
            headAngNorm.(wormcats{wormcatCtr}) = cell(numFiles,1);
            headAngSpeed.(wormcats{wormcatCtr}) = cell(numFiles,1);
            frameRunLengths.(wormcats{wormcatCtr}) = cell(numFiles,1);
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
            % apply phase restriction
            [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr);
            phaseFrameLogInd = trajData.frame_number <= lastFrame & trajData.frame_number >= firstFrame;
            trajData.filtered(~phaseFrameLogInd) = false;
            % find worms that have just left a cluster vs lone worms
            [leaveClusterLogInd, loneWormLogInd,~,~] = findWormCategory(filename,inClusterNeighbourNum,minNeighbrDist,postExitDuration);
            
            
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
                wormLogInd = trajData.worm_index_manual==uniqueWorms(wormCtr);
                wormtraj_xcoords = worm_xcoords(wormLogInd);
                wormtraj_ycoords = worm_ycoords(wormLogInd);
                fullTraj = figure; hold on
                if ~colorSpeed
                    plot(wormtraj_xcoords,wormtraj_ycoords,'k')
                else
                    nullheights = zeros(size(wormtraj_xcoords));
                    colors = midbodySpeedSigned(wormLogInd);
                    surface([wormtraj_xcoords,wormtraj_xcoords],...
                        [wormtraj_ycoords,wormtraj_ycoords],[nullheights,nullheights],[colors,colors],...
                        'facecol','no',...
                        'edgecol','interp',...
                        'linew',2);
                    colormap('jet')
                    caxis([0 10])
                end

                % now plot individually categorised traj on top of the full traj
                wormLogInd = trajData.worm_index_manual==uniqueWorms(wormCtr) & trajData.filtered;
                for wormcatCtr = 1:length(wormcats)
                    wormCatLogInd = wormLogInd & eval([wormcats{wormcatCtr} 'LogInd']);
                    if ~colorDirection
                        wormtraj_xcoords = worm_xcoords(wormCatLogInd);
                        wormtraj_ycoords = worm_ycoords(wormCatLogInd);
                        set(0,'CurrentFigure',fullTraj)
                        if strcmp(wormcats{wormcatCtr},'leaveCluster')
                            plot(wormtraj_xcoords,wormtraj_ycoords,'r.','MarkerSize',10)
                        elseif strcmp(wormcats{wormcatCtr},'loneWorm')
                            plot(wormtraj_xcoords,wormtraj_ycoords,'b.','MarkerSize',10)
                        end
                    else
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
                if ~colorDirection
                    figurename = (['figures/turns/redManualTrajFull/redManualTrajFull_' phase '_' filename(end-31:end-18) '_worm' num2str(uniqueWorms(wormCtr))]);
                else
                    figurename = (['figures/turns/redManualTrajFull/redManualTrajFullDir_' phase '_' filename(end-31:end-18) '_worm' num2str(uniqueWorms(wormCtr))]);
                end
                exportfig(fullTraj,[figurename '.eps'],exportOptions)
                system(['epstopdf ' figurename '.eps']);
                system(['rm ' figurename '.eps']);
            end
        end
    end
end