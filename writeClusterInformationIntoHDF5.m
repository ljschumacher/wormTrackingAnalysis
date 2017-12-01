% calculate distances to nearest neighbr and number of close neighbrs and
% write these into the hdf5 files

clear
close all

pixelsize = 100/19.5; % 100 microns are 19.5 pixels

strains = {'N2','npr1'};
wormnums = {'40'};
intensityThresholds_g = containers.Map({'40','HD','1W'},{60, 40, 100});
maxBlobSize_r = 2.5e5;
maxBlobSize_g = 1e4;
minSkelLength = 850;
maxSkelLength = 1500;
clusterCutOff = 500;

for strainCtr = 1:length(strains)
    for wormnum = wormnums
        %% load data
        filenames_r = importdata(['datalists/' strains{strainCtr} '_' wormnum{1} '_r_list.txt']);
        filenames_g = importdata(['datalists/' strains{strainCtr} '_' wormnum{1} '_g_list.txt']);
        numFiles = length(filenames_g);
        for fileCtr = 1:numFiles % can be parfor
            filename_r = filenames_r{fileCtr};
            trajData_r = h5read(filename_r,'/trajectories_data');
            blobFeats_r = h5read(filename_r,'/blob_features');
            skelData_r = h5read(filename_r,'/skeleton');
            filename_g = filenames_g{fileCtr};
            trajData_g = h5read(filename_g,'/trajectories_data');
            % filter green channel data by blob size and intensity
            blobFeats_g = h5read(filename_g,'/blob_features');
            trajData_g.filtered = filterIntensityAndSize(blobFeats_g,pixelsize,...
                intensityThresholds_g(wormnum{1}),maxBlobSize_g);
            % filter red channel data by blob size and intensity
            if contains(filename_r,'55')||contains(filename_r,'54')
                intensityThreshold_r = 80;
            else
                intensityThreshold_r = 40;
            end
            trajData_r.filtered = filterIntensityAndSize(blobFeats_r,pixelsize,...
                intensityThreshold_r,maxBlobSize_r);
            % filter by skeleton length
            trajData_r.filtered = trajData_r.filtered&logical(trajData_r.is_good_skel)&...
                filterSkelLength(skelData_r,pixelsize,minSkelLength,maxSkelLength);
            %% calculate stats - red channel files
            try
                neighbr_distances = h5read(filename_r,'/neighbr_distances');
                min_neighbr_dist = h5read(filename_r,'/min_neighbr_dist');
                num_close_neighbrs = h5read(filename_r,'/num_close_neighbrs');
            catch
                disp(['Calculating cluster status from ' filename_r ])
                [ neighbr_distances, min_neighbr_dist, num_close_neighbrs ] ...
                    = calculateNeighbrDistance(trajData_r,trajData_g,pixelsize,clusterCutOff);
                % check lengths
                assert(size(neighbr_distances,1)==length(trajData_r.frame_number))
                assert(size(neighbr_distances,2)==10)
                assert(length(min_neighbr_dist)==length(trajData_r.frame_number))
                assert(length(num_close_neighbrs)==length(trajData_r.frame_number))
                % write stats to hdf5-file
                h5create(filename_r,'/neighbr_distances',size(neighbr_distances),...
                    'Datatype','single')
                h5write(filename_r,'/neighbr_distances',single(neighbr_distances))
                h5create(filename_r,'/min_neighbr_dist',...
                    size(min_neighbr_dist),'Datatype','single')
                h5write(filename_r,'/min_neighbr_dist',...
                    single(min_neighbr_dist))
                h5create(filename_r,'/num_close_neighbrs',...
                    size(num_close_neighbrs),'Datatype','uint16')
                h5write(filename_r,'/num_close_neighbrs',...
                    uint16(num_close_neighbrs))
            end
            % calculate stats - green channel files
            try
                neighbr_distances = h5read(filename_g,'/neighbr_distances');
                min_neighbr_dist = h5read(filename_g,'/min_neighbr_dist');
                num_close_neighbrs = h5read(filename_g,'/num_close_neighbrs');
            catch
                disp(['Calculating cluster status from ' filename_g ])
                [ neighbr_distances, min_neighbr_dist, num_close_neighbrs ] ...
                    = calculateNeighbrDistance(trajData_g,trajData_r,pixelsize,clusterCutOff);
                % check lengths
                assert(size(neighbr_distances,1)==length(trajData_g.frame_number))
                assert(size(neighbr_distances,2)==10)
                assert(length(min_neighbr_dist)==length(trajData_g.frame_number))
                assert(length(num_close_neighbrs)==length(trajData_g.frame_number))
                % write stats to hdf5-file
                h5create(filename_g,'/neighbr_distances',size(neighbr_distances),...
                    'Datatype','single')
                h5write(filename_g,'/neighbr_distances',single(neighbr_distances))
                h5create(filename_g,'/min_neighbr_dist',...
                    size(min_neighbr_dist),'Datatype','single')
                h5write(filename_g,'/min_neighbr_dist',...
                    single(min_neighbr_dist))
                h5create(filename_g,'/num_close_neighbrs',...
                    size(num_close_neighbrs),'Datatype','uint16')
                h5write(filename_g,'/num_close_neighbrs',...
                    uint16(num_close_neighbrs))
            end
        end
    end
end

%% first data set - green channel files only, one more strain, and slightly different filters

strains = {'npr1','HA','N2'};
wormnums = {'40','HD'};
intensityThresholds_g('40') = 50;

for strainCtr = 1:length(strains)
    for wormnum = wormnums
        %% load data
        filenames_g = importdata(['datalists/' strains{strainCtr} '_' wormnum{1} '_list.txt']);
        numFiles = length(filenames_g);
        for fileCtr = 1:numFiles % can be parfor
            filename_g = filenames_g{fileCtr};
            trajData_g = h5read(filename_g,'/trajectories_data');
            % filter green channel data by blob size and intensity
            blobFeats_g = h5read(filename_g,'/blob_features');
            trajData_g.filtered = filterIntensityAndSize(blobFeats_g,pixelsize,...
                intensityThresholds_g(wormnum{1}),maxBlobSize_g);
            %% calculate stats - green channel files
            try
                neighbr_distances = h5read(filename_g,'/neighbr_distances');
                min_neighbr_dist = h5read(filename_g,'/min_neighbr_dist');
                num_close_neighbrs = h5read(filename_g,'/num_close_neighbrs');
            catch
                disp(['Calculating cluster status from ' filename_g ])
                [ neighbr_distances, min_neighbr_dist, num_close_neighbrs ] ...
                    = calculateNeighbrDistance(trajData_g,[],pixelsize,clusterCutOff);
                % check lengths
                assert(size(neighbr_distances,1)==length(trajData_g.frame_number))
                assert(size(neighbr_distances,2)==10)
                assert(length(min_neighbr_dist)==length(trajData_g.frame_number))
                assert(length(num_close_neighbrs)==length(trajData_g.frame_number))
                % write stats to hdf5-file
                h5create(filename_g,'/neighbr_distances',size(neighbr_distances),...
                    'Datatype','single')
                h5write(filename_g,'/neighbr_distances',single(neighbr_distances))
                h5create(filename_g,'/min_neighbr_dist',...
                    size(min_neighbr_dist),'Datatype','single')
                h5write(filename_g,'/min_neighbr_dist',...
                    single(min_neighbr_dist))
                h5create(filename_g,'/num_close_neighbrs',...
                    size(num_close_neighbrs),'Datatype','uint16')
                h5write(filename_g,'/num_close_neighbrs',...
                    uint16(num_close_neighbrs))
            end
        end
    end
end