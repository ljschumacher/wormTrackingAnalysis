% calculate distances to nearest neighbr and number of close neighbrs and
% write these into the hdf5 files

clear
close all

pixelsize = 100/19.5; % 100 microns are 19.5 pixels

strains = {'N2'}%,'N2'};
wormnums = {'HD'}%,'HD'};
intensityThresholds_g = [100, 60, 40];
maxBlobSize = 2.5e5;
maxBlobSize_g = 1e4;
minSkelLength = 850;
maxSkelLength = 1500;

for strainCtr = 1:length(strains)
    for numCtr = 1:length(wormnums)
        wormnum = wormnums{numCtr};
        %% load data
        filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_r_list.txt']);
        filenames_g = importdata(['datalists/' strains{strainCtr} '_' wormnum '_g_list.txt']);
        numFiles = length(filenames);
        for fileCtr = 1:numFiles % can be parfor
            filename = filenames{fileCtr};
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            skelData = h5read(filename,'/skeleton');
            filename_g = filenames_g{fileCtr};
            trajData_g = h5read(filename_g,'/trajectories_data');
            % filter green channel data by blob size and intensity
            blobFeats_g = h5read(filename_g,'/blob_features');
            trajData_g.filtered = filterIntensityAndSize(blobFeats_g,pixelsize,...
                intensityThresholds_g(numCtr),maxBlobSize_g);
            % filter red channel data by blob size and intensity
            if contains(filename,'55')||contains(filename,'54')
                intensityThreshold = 80;
            else
                intensityThreshold = 40;
            end
            trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                intensityThreshold,maxBlobSize);
            % filter by skeleton length
            trajData.filtered = trajData.filtered&logical(trajData.is_good_skel)&...
                filterSkelLength(skelData,pixelsize,minSkelLength,maxSkelLength);
            %% calculate stats - red files
            try
                neighbr_distances = h5read(filename,'/neighbr_distances');
                min_neighbr_dist = h5read(filename,'/min_neighbr_dist');
                num_close_neighbrs = h5read(filename,'/num_close_neighbrs');
            catch
                disp(['Calculating cluster status from ' filename ])
                [ neighbr_distances, min_neighbr_dist, num_close_neighbrs ] ...
                    = calculateClusterStatus(trajData,trajData_g,pixelsize,500);
                % check lengths
                assert(size(neighbr_distances,1)==length(trajData.frame_number))
                assert(size(neighbr_distances,2)==10)
                assert(length(min_neighbr_dist)==length(trajData.frame_number))
                assert(length(num_close_neighbrs)==length(trajData.frame_number))
                % write stats to hdf5-file
                h5create(filename,'/neighbr_distances',size(neighbr_distances),...
                    'Datatype','single')
                h5write(filename,'/neighbr_distances',single(neighbr_distances))
                h5create(filename,'/min_neighbr_dist',...
                    size(min_neighbr_dist),'Datatype','single')
                h5write(filename,'/min_neighbr_dist',...
                    single(min_neighbr_dist))
                h5create(filename,'/num_close_neighbrs',...
                    size(num_close_neighbrs),'Datatype','uint16')
                h5write(filename,'/num_close_neighbrs',...
                    uint16(num_close_neighbrs))
            end
        end
    end
end