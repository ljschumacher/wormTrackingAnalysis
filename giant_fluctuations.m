close all; clear; clc;

filename = 'recording34.2g_X1_skeletons.hdf5'
% Read in the data first
track_data = h5read(filename, '/trajectories_data');

% Then split it to obtain the important information
x_data = track_data.coord_x;
y_data = track_data.coord_y;
frames = track_data.frame_number;

% How many frames are we interested in?
final_t = max(frames);
final_t = 1000;

% Want to be able to draw a box
% Specify box position as a vector of 4 elements [x y w h]
box_w = 1000
box_h = 600
box_pos = [500 700 box_w box_h];

box_pos(1) = round(rand(1)*(peak2peak(x_data)-box_pos(3)))+min(x_data);
box_pos(2) = round(rand(1)*(peak2peak(y_data)-box_pos(4)))+min(y_data);


worm_counts = []
figure;
% For each frame of interest...
 for t = 0:final_t
     subplot(1,2,1)
     hold on
     rectangle('Position', box_pos)
     axis equal
     xlim([min(x_data) max(x_data)])
     ylim([min(y_data) max(y_data)])

    % Find all entries with data for this frame
    t_indexes = find(frames==t);
    % and count how many of these entries there are
    num_worms = length(t_indexes);

    % Initialise empty matrices to store location coordinates
    coords = zeros(num_worms,2); 

    % Pick out the xy coordinates for these worms at this time
    for i = 1:num_worms
        coords(i,1) = x_data(t_indexes(i));
        coords(i,2) = y_data(t_indexes(i));
    end
    A = coords(:,1);
    B = coords(:,2);
    op =find(A>box_pos(1) & A<box_pos(1)+box_pos(3));
    
    counter = 0;
    for i = 1:length(op)
        if B(op(i)) > box_pos(2) & B(op(i)) < box_pos(2)+box_pos(4)
            scatter(A(op(i)), B(op(i)));
            counter = counter+1;
        end
    end
    worm_counts(end+1) = counter;
    
    hold off
    
    subplot(1,2,2)
    plot(0:t, worm_counts)
    xlabel('Frame')
    ylabel('Worms in box')
    pause(0.1)
    
 end