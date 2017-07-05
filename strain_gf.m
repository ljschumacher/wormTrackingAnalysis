%% Script for computing and comparing giant fluctuations for strains
% Takes a cell array of .txt files as the 'in_list', which contain the
% locations of the appropriate .hdf5 files to read in.

% Finding the appropriate videos for each strain N2, HA and npr-1
close all; clear; clc;

% Supply the list files, which contain the directory locations of the
% videos for the titular strains
in_list = {'N2_40_list.txt', 'HA_40_list.txt', 'npr1_40_list.txt'}

exp_store = {}
models = {}

% For each one of these list files provided i.e. for each group of videos..
for strain = 1:length(in_list)
    
    % --------------------- Finding the right files --------------------- %
    
    strain_name = in_list{strain}
    % Initialise empty cell array for storing the video locations
    file_locations = {}
    
    % Open the .txt file with the appropriate name
    my_list_file = fopen(strain_name)
    
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
    
    num_movies = length(file_locations)
    
    batch_x = [];
    batch_y = [];
    exponents = [];
    
    % For each movie file identified in this list...
    for m = 1:num_movies
        
        m
        movie_x = [];
        movie_y = [];
        
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
        % Wish to ignore the first and last 10% of the video
        cut_out = [round(final_t*0.1) final_t-round(final_t*0.1)];
        
        % How many frames within this window do we wish to sample
        to_sample = 150;
        num_reps = 20;
        
        % Set up intermediate box scale factors to iterate through
        scale = [0.1:0.1:1];
        

        % ------------------------ Calculating gf ----------------------- %
        
        % For each box scale to be trialled
        for size = 1:length(scale)
            
            % Want to be able to draw a box
            box_w = 1000;
            box_h = 1000;
            % Specify box position as a vector of 4 elements [x y w h]
            box_pos = [500 500 box_w box_h];
            % Scale box size to investigate effect of changing this on the gf
            box_pos = box_pos.*scale(size);
            
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
                data_store(1,n,size)  = mean(worm_counts);
                movie_x(end+1) = mean(worm_counts);
                batch_x(end+1) = mean(worm_counts);
                
                data_store(2,n,size)  = var(worm_counts);
                movie_y(end+1) = var(worm_counts);
                batch_y(end+1) = var(worm_counts);
            end
            
        end
        
        % First figure plotting ln(N) against ln(var) for each video
        figure(1+2*(strain-1))
        
        subplot(num_movies,1, m)
        hold on
        for sca = 1:length(scale)
            
            x = log(data_store(1,:,sca));
            y = log(data_store(2,:,sca));
            x = x(isfinite(x));
            y = y(isfinite(y));
            
            scatter(x, y);
        end
        
        xlabel('ln(N)')
        ylabel('ln(var)')
        
        title(filename, 'Interpreter', 'none')
        hold off
        pause(0.01)
        
        % Also derive the linear regression through these logged points
        x = log(movie_x);
        y = log(movie_y);
        x = x(isfinite(x));
        y = y(isfinite(y));
        
        screen_indices = x > -2;
        x = x(screen_indices);
        y = y(screen_indices);
        
        % Use the mldivide operator to solve for a
        a = x'\y';
        % Store this power in the exponent array
        exponents(end+1) = a;
        
        
        
    end
    
    % ------------------------ Plotting results ------------------------- %
    % Second figure: plots all [N,var] data for entire strain
    figure;
    
    x = batch_x;
    y = batch_y;    
    
    x = x(isfinite(x));
    y = y(isfinite(y));
    
    screen_indices = x > -2;
    x = x(screen_indices);
    y = y(screen_indices);
    
    a = x'\y';
    y_calc = a.*x;
    
    % Plot this line in addition to the refline y = x
    hold on  
    scatter(x, y)
    plot(x,y_calc, 'r') 
    
    
    % Also wish to improve the fit by including a y intercept, b.
    % Can calculate b by padding x with a column of ones and using \
    X1 = [ones(length(x),1) x'];
    b = X1\y';
    
    % Plot this y-adjusted line in chosen color
    plot(x, y_calc+b(1), 'k');
    hold off
    
    % Add plot details and a reference line y = x
    xlabel('N')
    ylabel('var')    
    
    %title(in_list{s}, 'Interpreter', 'none')
    refline(1,0)
    
    % Store the current eponent for the video, for analysis later
    exp_store{end+1} = exponents;
    
    
    % ---------------------- Fitting a linear model --------------------- %
    md1 = fitlm(x,y,'linear');
    models{end+1} = md1;
    figure; plot(md1)
    xlabel('ln(N)')
    ylabel('ln(var)')
    %title(in_list{s}, 'Interpreter', 'none')
    
    
    % Construct figure where regressions for all strains are plotted
    % together, without scatters of the underlying datapoints
    figure(19)
    hold on
    plot(x,y_calc+b(1))
    legend(in_list, 'Interpreter', 'none')
    xlabel('ln(N)')
    ylabel('ln(var)')
    
    if strain == 1
       x1 = x;
       y1 = y;
    elseif strain == 2
        x2 = x;
        y2 = y;
    else
        x3 = x;
        y3 = y;
    end
    
end

% Consider the mean exponents for the regressions calculated
mean_exp = [mean(exp_store{1}), mean(exp_store{2}),mean(exp_store{3})]
hold off


%% Optional string manipulation of data to convert to R appropriate format
% x = x3; 
% y = y3;
% 
% my_string = 'data: ['
% for i = 1:length(x)
%     my_string = strcat(my_string, '[',sprintf('%1.1f',x(i)),', ', sprintf('%1.1f',y(i)),'], ');
% end
% my_string(end) = ''
% my_string = strcat(my_string,']')
% 
% x1try = x1;
% x2try = x2;
% x3try = x3;
% y1try = y1;
% y2try = y2;
% y3try = y3;
% 
% for i = 1:length(x1try)
%     x1try(i) = str2num(sprintf('%1.1f',x1try(i)));
% end
% for i = 1:length(x2try)
%     x2try(i) = sprintf('%1.1f',x2try(i));
% end
% for i = 1:length(x3try)
%     x3try(i) = sprintf('%1.1f',x3try(i));
% end
% for i = 1:length(y1try)
%     y1try(i) = sprintf('%1.1f',y1try(i));
% end
% for i = 1:length(y2try)
%     y2try(i) = sprintf('%1.1f',y2try(i));
% end
% for i = 1:length(y3try)
%     y3try(i) = sprintf('%1.1f',y3try(i));
% end