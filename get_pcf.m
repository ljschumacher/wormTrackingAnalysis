function [data_store] = get_pcf(filename, num_slices, slice_size)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to read in a .hdf5 file of choice and produce summary plots of 
% the pair correlation funtion g(r) for sampled frames at different points
% in the input video.

% Requires
% shadedErrorBar.m - Matlab function script to plot aesthetic shaded error
%                    bars around a mean. Must be present in current working
%                    directory. Available at https://uk.mathworks.com/matla
%                    bcentral/fileexchange/26311-raacampbell-shadederrorbar

% Inputs
% filename - string of the filename to be opened, contained in single 
%            quotes and complete with .hdf5 extension i.e. 
%            'recording34.2g_X1_skeletons.hdf5'.
% num_slices - integer for the number of smaller slices to be taken through
%              the video for more efficient analysis. Must be >= 3. These 
%              slices are taken at the very start and end of the input 
%              video, and then equally distributed inbetween.
% slice_size - non-zero, positive integer value defining the width of these
%              slices in frames.

% Outputs
% Produce three successive figures, presenting and comparing g(r). Also
% outputs a n by m matrix of pair correlation function values for all n
% frames of the video and for the m distance bins.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc;

max_gr_store = {};
r_pos_store = {};

for movie = 1:length(filename)
    % Read in the data first
    track_data = h5read(filename{movie}, '/trajectories_data');

    % Then split it to obtain the important information
    x_data = track_data.coord_x;
    y_data = track_data.coord_y;
    frames = track_data.frame_number;

    % How many frames are we interested in?
    final_t = max(frames);

    % Initialise array for storing the mean distances
    mean_dist = [];

    % For each frame of interest...
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
        bin_width = 0.5;
        bins = 1:bin_width:max(peak2peak(x_data)*p2m, peak2peak(y_data)*p2m); 
        %bins = 1:bin_width:max(peak2peak(pair_dist), peak2peak(pair_dist)); 
        
        % Get the histogram counts of the pair_dist data using the bins
        gr1 = histcounts(pair_dist,bins,'Normalization','pdf');
        gr1 = histcounts(pair_dist,bins);%,'Normalization','pdf');

        % Radial distribution function
        % Normalization step
        R = max(bins);
        gr2 = gr1.*R^2./(2*bins(2:end)*bin_width*((num_worms^2)-num_worms));

        % Store all data in a cell array for picking apart later
        if t == 0
            data_store = zeros(final_t+1, length(gr2));
        end
        data_store(t+1,:) = gr2;

        % Progress indicator
        Progress_reading_frames = single(t)/single(final_t)
    end


    % Consider only smaller parts of the video by sampling individual chunks
    % First initialise an empty cell array for storing legend strings.
    line_names = {};

    % Within the window constructed by chunks, consider every nth frame
    n = 10;

    % Obtain a matrix of the appropriate indexes for starting and ending
    % each chunk. These are uniformly split along the video, with a chunk at
    % the very start and end too
    chunks = [0 slice_size];
    if num_slices >= 3
        for i = 2:num_slices-1
            chunks(i,1) = round((i-1)*((final_t-(2*slice_size))/(num_slices-1))+(slice_size/2));
            chunks(i,2) = round((i-1)*((final_t-(2*slice_size))/(num_slices-1))+((3*slice_size)/2));
        end
    end
    chunks(num_slices,:) = [final_t-slice_size, final_t];

    % Set up a pallete for colouring the lines later
    colormap cool  
    pallete = colormap(1);
    col_size = size(pallete);

    % Now access the appropriate data using these chunk indices
    figure;


    % Consider how g(rmax) varies over time
    % This would be useful for comparing different movies to eachother

    % Initialise empty array to store values in
    max_gr = zeros(1,final_t+1);
    r_pos = zeros(1, final_t+1);

    % For each frame in the video
    for i = 1:final_t+1
        [max_gr(i), r_pos(i)] = max(data_store(i,:));
    end

    subplot(1,2,1)
    hold on
    plot((0:final_t),max_gr, 'Color', [0.4,0.4,0.4])
    xlabel('Frames')
    ylabel('g(r_{max})')

    % For each vertical divisor to be plotted
    for i = 1:num_slices
        % Get the appropriate colour index
        col_index = round(((col_size(1)-1)*mean(chunks(i,:)))/final_t);

        % Draw the vertical lines
        line([chunks(i,1),chunks(i,1)],ylim, 'Color', pallete(col_index+1,:))
        line([chunks(i,2),chunks(i,2)],ylim, 'Color', pallete(col_index+1,:))
    end
    hold off

    % Using the shadedErrorBar.m script to plot the variability in g(r) for 
    % each chunk of the video
    subplot(1,2,2)

    for i = 1:num_slices

        % Pick out the right appropriate g(r) data and assign to x and y
        indexes = chunks(i,:);
        y=data_store(indexes(1)+1:indexes(2)+1,:); x=bins(1:end-1);

        % Colour the line appropriately
        col_index = round(((col_size(1)-1)*mean(chunks(i,:)))/final_t);

        % Plot the line and standard deviation
        shadedErrorBar(x,y,{@mean,@std}, {'Color', pallete(col_index+1,:)}, 1);

        % Get the appropriate legend name
        line_names{i} = ['Frames ', num2str(chunks(i,1)), '-', num2str(chunks(i,2))];

        % Overlay transparent lines by calling hold on
        hold on
    end

    % Provide appropriate labels for axes
    xlabel('r [millimetres]')
    ylabel('g(r)')

    xlim([0.5 R-bin_width])

    % Finally add labelled colorbar and hold off to publish the plot
    set(gca, 'CLim', [0, single(final_t)])
    c = colorbar();
    c.Label.String = 'Frames from video';
    colormap cool  
    %legend(line_names)
    suptitle(strrep(filename{movie},'_','\_'))
    hold off
    
    
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
max_gr_store{end+1} = max_gr(:)
r_pos_store{end+1} = r_pos(:)

end

% Perform the appropriate plotting, having considered all input movies
figure;
subplot(2,1,1)
hold on
for movie = 1:length(filename)
    
    plot((0:final_t),smooth(max_gr_store{movie}',50))
    xlabel('Frames')
    ylabel('g(r_{max})')

end
legend(filename, 'Interpreter', 'none')
hold off 

subplot(2,1,2)
hold on
for movie = 1:length(filename)
    
    plot((0:final_t),smooth(r_pos_store{movie}',50))
    xlabel('Frames')
    ylabel('r_{max}') 
    
end
hold off 

end
