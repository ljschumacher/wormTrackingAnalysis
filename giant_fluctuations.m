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
% final_t = 10000;

% Want to be able to draw a box
box_w = 2000
box_h = 2000
% Get default line colors used by matlab and store as matrix
line_cols =  [0  0.4470 0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;]

% Decide how many repeat boes to draw at each of the scale factors
repeat_num = 7;
scale = [0.2:0.2:1]

% Initialise 3D matrix to store count data
worm_counts = zeros(repeat_num, final_t+1, length(scale));
figure;

% For each intermediate box size

for s = 1:length(scale)
    % Specify box position as a vector of 4 elements [x y w h]
    box_pos = [500 700 box_w box_h];
    % Scale box size to investigate effect of changing this on the gf
    box_pos = box_pos.*scale(s);
    s
    % For each repeat to be made
    for n = 1:repeat_num
 
        % Spawn a new box in the observed space
        box_pos(1) = round(rand(1)*(peak2peak(x_data)-box_pos(3)))+min(x_data);
        box_pos(2) = round(rand(1)*(peak2peak(y_data)-box_pos(4)))+min(y_data);

         % For each frame of interest...
         for t = 0:final_t

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
            
            % Get the indices for the coordinates in the correct x range
            % for the current box
            op =find(A > box_pos(1) & A < box_pos(1) + box_pos(3));
            
            % Count how many of these suitable x coordinates also have
            % appropriate y coordinate pairs
            counter = 0;
            for i = 1:length(op)
                if B(op(i)) > box_pos(2) & B(op(i)) < box_pos(2)+box_pos(4)
                    counter = counter+1;
                end
            end
            % Adjust the count for the number of worms in the box suitably
            worm_counts(n,t+1,s) = counter;

         end
    
    subplot(2,length(scale),s+length(scale))
    hold on
    rectangle('Position', box_pos, 'EdgeColor', line_cols(n,:))
    axis equal
    xlim([min(x_data) max(x_data)])
    ylim([min(y_data) max(y_data)])
    hold off     
    end


    subplot(2,length(scale),s)
    hold on
    for n = 1:repeat_num
        plot(0:t, smooth(worm_counts(n,:,s),30))
    end
    hold off
    xlabel('Frame')
    ylabel('Worms in box')
end

% Perform some work examing the variance from this frequency data
vars = zeros(length(scale), final_t+1);
mean_vars = []

for i = 1:length(scale)

    for t = 1:final_t+1
        vars(i,t) = var(worm_counts(:,t,i));
    end
    mean_vars(end+1) = mean(vars(i,:))
end

figure;
hold on
plot(scale, mean_vars)
hold off
xlabel('Box scaling factor')
ylabel('Mean variation in worm frequencies')
