% calculate various single worm statistics for different numbers of worms
% on plate

% issues / todo:
% - should distances be calculated only to other moving worms, or to any
% object (ie also worms that won't appear in the next frame)?

clear
close all

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

pixelsize = 100/19.5; % 100 microns are 19.5 pixels

strains = {'npr1','HA','N2'};
wormnums = {'40','HD','1W'};
midbodyIndcs = 19:33;


for strainCtr = 1:length(strains)
    revFreqFig = figure; hold on
    revDurFig = figure; hold on
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        %% load data
        filenames = importdata([strains{strainCtr} '_' wormnum '_r_list.txt']);
        numFiles = length(filenames);
        midbodySpeeds = cell(numFiles,1);
        reversalfreq_lone = NaN(numFiles,1);
        reversaldurations_lone = cell(numFiles,1);
        for fileCtr = 1:numFiles % can be parfor
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            % %             featData = h5read(strrep(filename,'skeletons','features'),'/features_timeseries');
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            maxNumFrames = numel(unique(trajData.frame_number));
            numFrames = maxNumFrames;
            framesAnalyzed = unique(trajData.frame_number); % randomly sample frames without replacement
            %% calculate stats
            if ~strcmp(wormnum,'1W')
                loneWorms = cell(numFrames,1);
                for frameCtr = 1:numFrames
                    frame = framesAnalyzed(frameCtr);
                    %%% need to update this for 2 channel case
                    [~, loneWorms{frameCtr}, ~] = ...
                        getWormClusterStatus(trajData, frame, pixelsize, 2500,500,2);
                end
                % pool data from all frames
                loneWorms = horzcat(loneWorms{:});
            else
                loneWorms = true(numFrames);
            end
            %% temporary: calculate overall midbody speeds, until we can link
            % trajectories and features from the tracker
            skelData = h5read(filename,'/skeleton');
            % centroids of midbody skeleton
            midbody_x = mean(squeeze(skelData(1,midbodyIndcs,:)))*pixelsize;
            midbody_y = mean(squeeze(skelData(2,midbodyIndcs,:)))*pixelsize;
            % change in centroid position over time (issue of worm shifts?)
            dmidbody_xdt = gradient(midbody_x)*frameRate;
            dmidbody_ydt = gradient(midbody_y)*frameRate;
            % midbody speed and velocity
            midbodySpeed = sqrt(dmidbody_xdt.^2 + dmidbody_ydt.^2)./gradient(double(trajData.frame_number))';
            midbodyVelocity = [dmidbody_xdt.^2; dmidbody_ydt]./gradient(double(trajData.frame_number))';
            % direction of segments pointing along midbody
            [~, dmidbody_yds] = gradient(squeeze(skelData(2,midbodyIndcs,:)),-1);
            [~, dmidbody_xds] = gradient(squeeze(skelData(1,midbodyIndcs,:)),-1);
            % sign speed based on relative orientation of velocity to midbody
            midbodySpeedSigned = sign(sum(midbodyVelocity.*[mean(dmidbody_xds); mean(dmidbody_yds)])).*midbodySpeed;
            % smooth speed to denoise
            midbodySpeedSigned = smooth(midbodySpeedSigned,3,'moving');
            % ignore first and last frames of each worm's track
            wormChangeIndcs = gradient(double(trajData.worm_index_joined))~=0;
            midbodySpeedSigned(wormChangeIndcs)=NaN;
            % find reversals in midbody speed
            [revStartInd, revDuration] = findReversals(...
                struct('midbody_speed',midbodySpeedSigned,'worm_index',trajData.worm_index_joined));
            loneReversals = ismember(revStartInd,find(loneWorms));
            Nrev_lone = nnz(loneReversals);
            T_lone = nnz(loneWorms)/frameRate;
            Trev_lone = nnz(midbodySpeedSigned(loneWorms)<0)/frameRate;
            reversalfreq_lone(fileCtr) = Nrev_lone./(T_lone - Trev_lone);
            reversaldurations_lone{fileCtr} = revDuration(loneReversals);
        end
        %% plot data
        reversaldurations_lone = vertcat(reversaldurations_lone{:});
        histogram(revDurFig.Children,reversaldurations_lone,0:(10*frameRate),'Normalization','pdf','EdgeColor','none');
      
        boxplot(revFreqFig.Children,reversalfreq_lone,numCtr)
    end
    %% format and export figures
    title(revFreqFig.Children,[strains{strainCtr} ', ' num2str(numSamples) ' samples'],'FontWeight','normal');
    set(revFreqFig,'PaperUnits','centimeters')
    %     xlabel(revFig.Children,'speed (\mum/s)')
    ylabel(revFreqFig.Children,'P')
    xlabel('frames')
    legend(wormnums)
    %     figurename = ['figures/singleWorm/' strains{strainCtr} '_speeddistributions'];
    exportfig(revFreqFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end