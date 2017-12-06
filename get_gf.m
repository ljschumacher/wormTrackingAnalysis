%% Script to visualise giant fluctuations for individual movies specified
% Takes a list of .hdf5 files and creates plots to show giant fluctuations
close all; clear; clc;

filename = {'recording34.2g_X1_skeletons.hdf5', 'recording34.5g_X1_skeletons.hdf5', 'recording34.8g_X1_skeletons.hdf5'}
num_movies = length(filename);

all_pairs_x = [];
all_pairs_y = [];

for movie = 1:num_movies
    % Initial filtering of data to eliminate larvae
    % First, set the established parameters used in Linus' GitHub.
    %intensityThresholds = containers.Map({'40', 'HD','1W'},{60,40,100});
    intensityThresholds = 50;
    maxBlobSize = 1e4;
    pixelsize = 100/19.5; % 100 microns are 19.5 pixels;
    
    % Read in the associated feature data
    blob_data = h5read(filename{movie}, '/blob_features');
    
    % Get filtering indices according to the parameters outlined
    filter = (blob_data.area*pixelsize^2 <= maxBlobSize) & ...
        (blob_data.intensity_mean>=intensityThresholds);
    
    
    % Read in the trajectory data
    track_data = h5read(filename{movie}, '/trajectories_data');
    
    % Then split it to obtain the important information, whilst
    % applying the filter obtained from \blob_features above
    x_data = track_data.coord_x(filter);
    y_data = track_data.coord_y(filter);
    frames = track_data.frame_number(filter);
    
    % Read in the data first
    track_data = h5read(filename{movie}, '/trajectories_data');
    
    % Then split it to obtain the important information
    x_data = track_data.coord_x;
    y_data = track_data.coord_y;
    frames = track_data.frame_number;
    
    % Set up intermediate box scale factors to iterate through. These
    % represent fractions of the larger field of view
    scale = [0.1:0.1:0.5];
    
    % Matrix to store the data
    num_reps = 20;
    data_store = zeros(2,num_reps, length(scale));
    
    % How many frames are we interested in?
    final_t = max(frames);
    % How many frames within this window do we wish to sample
    to_sample = 100;
    
    % Wish to ignore the first and last 10% of the video
    cut_out = [round(final_t*0.1) final_t-round(final_t*0.1)];
    
    
    % For each box scale to be trialled
    for s = 1:length(scale)
        
        % Want to be able to draw a box
        box_w = 1000;
        box_h = 1000;
        % Specify box position as a vector of 4 elements [x y w h]
        box_pos = [500 500 box_w box_h];
        % Scale box size to investigate effect of changing this on the gf
        box_pos = box_pos.*scale(s);

        for n = 1:num_reps
            % Spawn a new box in the observed space
            box_pos(1) = round(rand(1)*(peak2peak(x_data)-box_pos(3)))+min(x_data);
            box_pos(2) = round(rand(1)*(peak2peak(y_data)-box_pos(4)))+min(y_data);            
            
            % Obtain these random frame indices within the defined limits
            rand_frames = randi(peak2peak(cut_out), to_sample, 1);
            rand_frames = rand_frames+double(cut_out(1));
            
            % Initialise matrix to store count values in later
            worm_counts = zeros(1, to_sample);
            
            for i = 1:length(rand_frames)
                t = rand_frames(i);
                % Find all entries with data for this frame
                t_indexes = find(frames==t);
                
                % and count how many of these entries there are
                num_worms = length(t_indexes);
                
                % Initialise empty matrices to store location coordinates
                coords = zeros(num_worms,2);
                
                % Pick out the xy coordinates for these worms at this time
                for j = 1:num_worms
                    coords(j,1) = x_data(t_indexes(j));
                    coords(j,2) = y_data(t_indexes(j));
                end
                A = coords(:,1);
                B = coords(:,2);
                
                % Get the indices for the coordinates in the correct x range
                % for the current box
                op =find(A > box_pos(1) & A < box_pos(1) + box_pos(3));
                
                % Count how many of these suitable x coordinates also have
                % appropriate y coordinate pairs
                counter = 0;
                for q = 1:length(op)
                    if B(op(q)) > box_pos(2) & B(op(q)) < box_pos(2)+box_pos(4)
                        counter = counter+1;
                    end
                end

                % Adjust the count for the number of worms in the box suitably
                worm_counts(1,i) = counter;
            end
            
            % Store the data in appropriate matrixes for plotting later
            data_store(1,n,s)  = mean(worm_counts);
            all_pairs_x(end+1) = mean(worm_counts);
            
            data_store(2,n,s)  = var(worm_counts);
            all_pairs_y(end+1) = var(worm_counts);
        end

    end
    
end

% Also plot points from videos together in separate fig & color accordingly
figure;
hold on
slice = length(all_pairs_x)/num_movies;
for movie = 0:num_movies-1
    scatter(all_pairs_x((1+slice*movie):(slice*(movie+1))),all_pairs_y((1+slice*movie):(slice*(movie+1))))
end
hold off

% Add plot details and a reference line y = x
xlabel('N')
ylabel('var')
refline(1,0)
legend(filename, 'Interpreter', 'none')


% Fitting linear relationships through the data

figure;
slice = length(all_pairs_x)/num_movies;

for movie = 0:num_movies-1
    x = all_pairs_x;
    y = all_pairs_y;
    
    x = x((1+slice*movie):(slice*(movie+1)));
    y = y((1+slice*movie):(slice*(movie+1)));
    
    x = x(isfinite(x));
    y = y(isfinite(y));

    b1 = x'\y';
    y_calc = b1.*x;
    
    subplot(3,1,movie+1)
    hold on
    scatter(x,y)
    plot(x,y_calc)
    refline(1,0)

    xlabel('N')
    ylabel('var')
    legend({'data', ['a = ' num2str(b1)], 'a = 1',})
    title(filename{movie+1}, 'Interpreter', 'none')
    hold off
end