%% Script to calculate and compare the pair correlation for different strains
close all; clear; clc;

% Supply the list files, which contain the directory locations of the
% videos for the titular strains
in_list = {'N2_40_list.txt', 'HA_40_list.txt', 'npr1_40_list.txt'};

exp_store = {};
models = {};
mean_store = {};
shade_store = {};

% For each one of these list files provided i.e. for each group of videos..
for strain = 1:length(in_list)
    
    max_gr_store = {};
    r_pos_store = {};
    % --------------------- Finding the right files --------------------- %
    
    strain_name = in_list{strain};
    % Initialise empty cell array for storing the video locations
    file_locations = {};
    
    % Open the .txt file with the appropriate name
    my_list_file = fopen(strain_name);
    
    % Read in the first line of the .txt file
    tline = fgetl(my_list_file);
    % For every line that still contains text...
    while ischar(tline)
        % Store this file location in the initialised array
        file_locations{end+1} = tline;
        % Read in the next line
        tline = fgetl(my_list_file);
    end
    % Close the list file, having read in all of the file locations
    fclose(my_list_file)
    
    num_movies = length(file_locations);
    
    % For each movie file identified in this list...
    for m = 1:num_movies       
        
        % Manipulate the stored string to get the filename
        full = file_locations{m};
        find(full=='/',1, 'last');
        short = full(find(full=='/',1, 'last')+1:end);
        
        % Store the short filename for opening the video
        filename = short;
        % Extract and store the strain id from the name i.e. '34.8'
        vid_id = short(10:13);
        
        % ----------------- Obtaining & Filtering data ------------------ %
        
        % Must obtain logical indices for filtering the tracking data based
        % on the intensity and size of potential worms. This is to
        % eliminate larvae from the files
        
        % First, set the established parameters used in Linus' GitHub.
        %intensityThresholds = containers.Map({'40', 'HD','1W'},{60,40,100});
        intensityThresholds = 50;
        maxBlobSize = 1e4;
        pixelsize = 100/19.5; % 100 microns are 19.5 pixels;
        
        % Read in the associated feature data
        blob_data = h5read(filename, '/blob_features');
        
        % Get filtering indices according to the parameters outlined
        filter = (blob_data.area*pixelsize^2 <= maxBlobSize) & ...
            (blob_data.intensity_mean>=intensityThresholds);
        
                
        % Read in the trajectory data
        track_data = h5read(filename, '/trajectories_data');
        
        % Then split it to obtain the important information, whilst
        % applying the filter obtained from \blob_features above
        x_data = track_data.coord_x(filter);
        y_data = track_data.coord_y(filter);
        frames = track_data.frame_number(filter);
        
        % How many frames are we interested in?
        final_t = max(frames);
        
        % Initialise array for storing the mean distances
        mean_dist = [];
        
        for t = 0:final_t

            % Find all entries with data for this frame
            t_indexes = find(frames==t);
            % and count how many of these entries there are
            num_worms = length(t_indexes);

            % Initialise empty matrices to store location coordinates
            coords = zeros(num_worms,2); 

            % Initialise empty matrix to store pairwise distances
            distances = [];

            % Pick out the xy coordinates for these worms at this time
            for i = 1:num_worms
                coords(i,1) = x_data(t_indexes(i));
                coords(i,2) = y_data(t_indexes(i));
            end

            % Obtain the pairwise distances with pdist
            pair_dist = pdist(coords);

            % Make the distances informative
            % There are 100 microns in a pixel
            p2m = 0.1/19.5;
            pair_dist = pair_dist.*p2m;
            mean_dist(end+1) = mean(pair_dist);

            % Create bins of a given width, to store the data in
            bin_width = 0.25;
            bins = 0:bin_width:max(peak2peak(x_data)*p2m, peak2peak(y_data)*p2m); 
            %bins = 1:bin_width:max(peak2peak(pair_dist), peak2peak(pair_dist)); 

            % Get the histogram counts of the pair_dist data using the bins
            gr1 = histcounts(pair_dist,bins,'Normalization','count');

            % Radial distribution function
            % Normalization step
            R = max(bins);
            gr2 = gr1.*R^2./(2*bins(2:end)*bin_width*((num_worms^2)-num_worms));

            % Store all data in a cell array for picking apart later
            if t == 0
                data_store = zeros(final_t+1, length(gr2));
            end
            data_store(t+1,:) = gr2;

        end
        
        
        % Consider how g(rmax) varies over time
        % This would be useful for comparing different movies to eachother
        % Initialise empty array to store values in
        max_gr = zeros(1,final_t+1);
        r_pos = zeros(1, final_t+1);
        
        % For each frame in the video
        for i = 1:final_t+1
            % Get the maximum value of g(r) in that row, and its index
            [max_gr(i), r_pos(i)] = max(data_store(i,:));
        end
        
        % Convert this index to the appropriate distance, using bin_width
        r_pos = (r_pos.*bin_width) - 0.5*bin_width;
        
        % Store these values in the initialised cell arrays for plotting
        max_gr_store{end+1} = max_gr(:);
        r_pos_store{end+1} = r_pos(:);
        
    end
    
    % Plotting section
    figure(10)
    hold on
    for movie = 1:num_movies
        ys = max_gr_store{movie}';
        if strain == 1
            plot((0:9999),smooth(ys(1:10000),50), 'Color','red')
        elseif strain == 2
            plot((0:9999),smooth(ys(1:10000),50), 'Color','green')
        else
            plot((0:9999),smooth(ys(1:10000),50), 'Color','blue')
        end
        xlabel('Frames')
        ylabel('g(r_{max})')
    end
    pause(0.01)
    
plotting_data = zeros(10000,num_movies);
for movie = 1:num_movies
    initial = max_gr_store{movie};
    plotting_data(:,movie) = initial(1:10000);
end

means = zeros(1,10000);
for i = 1:length(means)
    means(i) = mean(plotting_data(i,:));
end   

shade_data = zeros(num_movies*100,100)

for q = 1:100
    trial = plotting_data'
    happy = trial(:,((q-1)*100)+1:((q-1)*100)+100)
    B = happy(:)
    shade_data(:,q) = B
end

mean_store{strain} = means
shade_store{strain} = shade_data

end


hold off

figure(20)
for i = 1:length(in_list)
    hold on
    plot(1:10000, smooth(mean_store{i},100))
end
hold off

xlabel('Frames')
ylabel('g(r_{max})')

figure(30)
shadedErrorBar((1:10000),plotting_data',{@mean,@std}, {'Color', 'red'})

figure(40)
for i = 1:length(in_list)
    hold on
    shade_data = shade_store{i};
    shade_data(isnan(shade_data))=0;
    
    if i == 1
        shadedErrorBar((1:100:10000),shade_data,{@mean,@std}, {'Color', 'red'}, 1)
    elseif i == 2
        shadedErrorBar((1:100:10000),shade_data,{@mean,@std}, {'Color', 'green'}, 1)
    else
        shadedErrorBar((1:100:10000),shade_data,{@mean,@std}, {'Color', 'blue'}, 1)
    end
end

hold off
xlabel('Frames')
ylabel('g(r_{max})')

figure(50)
shadedErrorBar((1:100:10000),shade_data,{@mean,@std}, {'Color', 'blue'})
