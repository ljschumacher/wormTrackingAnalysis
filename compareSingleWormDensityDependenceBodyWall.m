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

strains = {'npr1','N2'};
wormnums = {'1W','40'};%,'HD'};
midbodyIndcs = 19:33;
plotColors = lines(length(wormnums));

for strainCtr = 1:length(strains)
    revFreqFig = figure; hold on
    revDurFig = figure; hold on
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        %% load data
        filenames = importdata([strains{strainCtr} '_' wormnum '_r_list.txt']);
        if ~strcmp(wormnum,'1W')
            filenames_g = importdata([strains{strainCtr} '_' wormnum '_g_list.txt']);
        else
            filenames_g = {};
        end
        numFiles = length(filenames);
        midbodySpeeds = cell(numFiles,1);
        reversalfreq_lone = NaN(numFiles,1);
        reversaldurations_lone = cell(numFiles,1);
        reversalfreq_incluster = NaN(numFiles,1);
        reversaldurations_incluster = cell(numFiles,1);
        for fileCtr = 1:numFiles % can be parfor?
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            if ~strcmp(wormnum,'1W')
                filename_g = filenames_g{fileCtr};
                trajData_g = h5read(filename_g,'/trajectories_data');
            end
            % %             featData = h5read(strrep(filename,'skeletons','features'),'/features_timeseries');
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            maxNumFrames = numel(unique(trajData.frame_number));
            numFrames = maxNumFrames;
            framesAnalyzed = unique(trajData.frame_number);
            %% calculate stats
            if ~strcmp(wormnum,'1W')
                loneWorms = cell(numFrames,1);
                inCluster = cell(numFrames,1);
                for frameCtr = 1:numFrames
                    frame = framesAnalyzed(frameCtr);
                    %%% need to update this for 2 channel case
                    [~, loneWorms_rr, ~] = ...
                        getWormClusterStatus(trajData, frame, pixelsize, 1300, 500, 3);
                    [inCluster{frameCtr}, loneWorms_rg, ~] = ...
                    getWormClusterStatus2channel(trajData, trajData_g, frame, pixelsize, 1900, 500, 3);
                    loneWorms{frameCtr} = loneWorms_rr&loneWorms_rg;
                end
                % pool data from all frames
                loneWorms = horzcat(loneWorms{:});
                inCluster = vertcat(inCluster{:});
            else
                loneWorms = true(size(trajData.frame_number));
                inCluster = false(size(trajData.frame_number));
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
            % ignore frames with bad skeletonization
            midbodySpeedSigned(trajData.is_good_skel~=1)=NaN;
            % find reversals in midbody speed
            [revStartInd, revDuration] = findReversals(...
                struct('midbody_speed',midbodySpeedSigned,'worm_index',trajData.worm_index_joined));
            loneReversals = ismember(revStartInd,find(loneWorms));
            inclusterReversals = ismember(revStartInd,find(inCluster));
            Nrev_lone = nnz(loneReversals);
            Nrev_incluster = nnz(inCluster);
            T_lone = nnz(loneWorms)/frameRate;
            T_incluster = nnz(inCluster)/frameRate;
            Trev_lone = nnz(midbodySpeedSigned(loneWorms)<0)/frameRate;
            Trev_incluster = nnz(midbodySpeedSigned(inCluster)<0)/frameRate;
            reversalfreq_lone(fileCtr) = Nrev_lone./(T_lone - Trev_lone);
            reversalfreq_incluster(fileCtr) = Nrev_incluster./(T_incluster - Trev_incluster);
            reversaldurations_lone{fileCtr} = revDuration(loneReversals)/frameRate;
            reversaldurations_incluster{fileCtr} = revDuration(inclusterReversals)/frameRate;
        end
        %% plot data
        reversaldurations_lone = vertcat(reversaldurations_lone{:});
        reversaldurations_incluster = vertcat(reversaldurations_incluster{:});
        histogram(revDurFig.Children,reversaldurations_lone,0:1/3:30,...
            'Normalization','pdf','DisplayStyle','stairs','EdgeColor',plotColors(numCtr,:));
        histogram(revDurFig.Children,reversaldurations_incluster,0:1/3:30,...
            'Normalization','pdf','EdgeColor','none','FaceColor',plotColors(numCtr,:));
        
        boxplot(revFreqFig.Children,reversalfreq_lone,'Positions',numCtr-1/4,...
            'Notch','off')
        boxplot(revFreqFig.Children,reversalfreq_incluster,'Positions',numCtr+1/4,...
            'Notch','off')
        revFreqFig.Children.XLim = [0 length(wormnums)+1];
    end
    %% format and export figures
    title(revFreqFig.Children,strains{strainCtr},'FontWeight','normal');
    set(revFreqFig,'PaperUnits','centimeters')
    revFreqFig.Children.XTick = 1:length(wormnums);
    revFreqFig.Children.XTickLabel = wormnums;
    revFreqFig.Children.XLabel.String = 'worm number';
    revFreqFig.Children.YLabel.String = 'reversals/time';
    revFreqFig.Children.YLim = [0 1];
    figurename = ['figures/singleWorm/' strains{strainCtr} '_reversals'];
    exportfig(revFreqFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
    %
    title(revDurFig.Children,strains{strainCtr},'FontWeight','normal');
    set(revDurFig,'PaperUnits','centimeters')
    xlabel(revDurFig.Children,'time (s)')
    ylabel(revDurFig.Children,'P')
    legend(revDurFig.Children,wormnums)
    figurename = ['figures/singleWorm/' strains{strainCtr} '_reversaldurations'];
    exportfig(revDurFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
end