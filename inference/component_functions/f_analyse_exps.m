function [exp_ss_array, in_list] = f_analyse_exps(in_list, using_second)
global num_statistics

% Construct arrays for storing the summary statistic outputs
exp_ss_array = cell(length(in_list), num_statistics);
strain_vars = zeros(length(in_list), num_statistics-1);

% For each one of these list files provided i.e. for each group of videos..
for strain = 1:length(in_list)
    
    % FINDING THE RIGHT FILES %
    strain_name = in_list{strain};
    exp_ss_array{strain,1} = strain_name;
    % Initialise empty cell array for storing the video locations
    exp_file_locations = {};
    
    % Open the .txt file with the appropriate name
    my_list_file = fopen(strain_name);
    
    % Read in the first line of the .txt file
    tline = fgetl(my_list_file);
    % For every line that still contains text...
    while ischar(tline)
        % Store this file location in the initialised array
        exp_file_locations{end+1} = tline;
        % Read in the next line
        tline = fgetl(my_list_file);
    end
    % Close the list file, having read in all of the file locations
    fclose(my_list_file)
    
    num_movies = length(exp_file_locations);
    
    exp_replicate_ss_array = cell(num_movies, num_statistics);
    % For each movie file identified in this list...
    for m = 1:num_movies
        
        % Manipulate the stored string to get the filename
        full = exp_file_locations{m};
        find(full=='/',1, 'last');
        short = full(find(full=='/',1, 'last')+1:end);
        
        % Store the short filename for opening the video
        filename = short;
        % Extract and store the strain id from the name i.e. '34.8'
        vid_id = short(10:13);
        
        % OBTAINING AND FILTERING DATA 
        
        % Must obtain logical indices for filtering the tracking data based
        % on the intensity and size of potential worms. This is to
        % eliminate larvae from the files
        
        % First, set the established parameters used in Linus' GitHub.
        if using_second == 0;
            intensityThresholds = 50;
        elseif using_second == 1;
            intensityThresholds = 60;
        end
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
        
        % Need to omit opening/specified frames to obtain stationary phase
        % Use xlsread to read in spreadsheet data containing end frames
        
        try
            [lastFrames, filenames, ~] = xlsread([strtok(strain_name,'.') '.xlsx'], 1, 'A1:B15','basic');
            lastFrame = lastFrames(m);
        catch
            if m == 1
                fprintf(['Could not find the list file ' strtok(strain_name,'.')...
                    '.xlsx \nWithout information on stationary phase, using maximum frame as endpoint instead.\n'])
            end
            last_Frame = max(frames);
        end
        
        % Omit the first 10% of frames and the associated worm coordinates
        firstFrame = double(round(max(track_data.frame_number)*0.1));
        stationary_filter = frames < lastFrame & frames > firstFrame;
        x_data = x_data(stationary_filter);
        y_data = y_data(stationary_filter);
        frames = frames(stationary_filter);
        
        in_data = {x_data, y_data, frames};
        
        % Then compute all chosen summary statistics, as with the simulated
        % data above.
        ss_results = f_compute_ss(in_data, 'experiment')
        
        for each_ss = 1:length(ss_results)
            exp_replicate_ss_array{m, each_ss+1} = ss_results{each_ss}
        end
                
    end
    
    % Combine SS from different experimental replicates, to give single
    % reference for each strain
    replicate_summary = cell(1,num_statistics-1);
    for elem = 2:num_statistics
        replicate_summary{elem-1} = mean(cell2mat(exp_replicate_ss_array(:,elem)));
        for mov = 1:num_movies
            strain_vars(strain,elem-1) = strain_vars(strain,elem-1) ... 
                + var(exp_replicate_ss_array{mov,elem});
        end
    end
    strain_vars(strain,:) = strain_vars(strain,:)/num_movies;
    exp_ss_array(strain, 2:end) = replicate_summary;
end
end