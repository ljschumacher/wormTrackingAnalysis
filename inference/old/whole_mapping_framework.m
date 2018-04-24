% A framework in which to create a mapping between simple and complex model
% space, in order to propose more informative priors for future complex
% simulation.

% Shortcut if sim_ss_array has already been made, saves recalculating for
% all 10000 simulations (3 hours)
sim_ss_defined = 1;
%% Compute the statistics of the simple simulations first

close all; clc;

% Firstly, read in the list of simulation files
simple_sim_file_list = 'simple_simulation_list_mapping.txt';
%sim_file_list = 'linus_sim_list.txt';
simple_sim_file_names = {};
my_list_file = fopen(simple_sim_file_list);

% Read in the first line of the .txt file
tline = fgetl(my_list_file);
% For every line that still contains text...
while ischar(tline)
    % Store this file location in the initialised array
    simple_sim_file_names{end+1} = tline;
    % Read in the next line
    tline = fgetl(my_list_file);
end
% Close the list file, having read in all of the file names
fclose(my_list_file)


% Initialise an array for storing the summary statistic outputs
num_statistics = 5;
simple_sim_ss_array = cell(length(simple_sim_file_names), num_statistics);

% For each simulation file in the list, compute the appropriate summary
% statistics using supplied functions
for sim = 1:length(simple_sim_file_names)
    computing_statistics_for_simple_simulations = sim/length(simple_sim_file_names)
    load(simple_sim_file_names{sim});
    simple_sim_ss_array{sim,1} = simple_sim_file_names{sim};
    % Read in the appropriate data
    data = xyarray;
    % Compute SS1: Distribution of the number of worms per cluster
    [SS1,SS2,SS3] = inf_clusterdistribution(data, 'simulation');
    simple_sim_ss_array{sim,2} = SS1;
    
    % SS2: Distribution of the branch heights from the same clustering
    simple_sim_ss_array{sim,3} = SS2;
    
    % SS3 : Percentage of worms in cluster
    simple_sim_ss_array{sim,4} = SS3;
    
    % SS4: gr array, showing distribution of neighbours with distance
    simple_sim_ss_array{sim,5} = inf_pcf(data, 'simulation');
    
end

%% Compute the statistics for the small number of complex simulations

% Firstly, read in the list of simulation files
complex_sim_file_list = 'complex_simulation_list_mapping.txt';
%sim_file_list = 'linus_sim_list.txt';
complex_sim_file_names = {};
my_list_file = fopen(complex_sim_file_list);

% Read in the first line of the .txt file
tline = fgetl(my_list_file);
% For every line that still contains text...
while ischar(tline)
    % Store this file location in the initialised array
    complex_sim_file_names{end+1} = tline;
    % Read in the next line
    tline = fgetl(my_list_file);
end
% Close the list file, having read in all of the file names
fclose(my_list_file)


% Initialise an array for storing the summary statistic outputs
num_statistics = 5;
complex_sim_ss_array = cell(length(complex_sim_file_names), num_statistics);

% For each simulation file in the list, compute the appropriate summary
% statistics using supplied functions
for sim = 1:length(complex_sim_file_names)
    progress_reading_complex_simulations = sim/length(complex_sim_file_names)
    load(complex_sim_file_names{sim});
    complex_sim_ss_array{sim,1} = complex_sim_file_names{sim};
    % Read in the appropriate data
    data = xyarray;
    % Compute SS1: Distribution of the number of worms per cluster
    [SS1,SS2,SS3] = inf_clusterdistribution(data, 'complexsim');
    complex_sim_ss_array{sim,2} = SS1;
    
    % SS2: Distribution of the branch heights from the same clustering
    complex_sim_ss_array{sim,3} = SS2;
    
    % SS3 : Percentage of worms in cluster
    complex_sim_ss_array{sim,4} = SS3;
    
    % SS4: gr array, showing distribution of neighbours with distance
    complex_sim_ss_array{sim,5} = inf_pcf(data, 'complexsim');
    
end









%* Then build the experimental reference with which to compare the simple
%  and complex simulations

% Supply the list files, which contain the directory locations of the
% videos for the titular strains

% First data set
in_list = {'N2_40_list.txt', 'HA_40_list.txt', 'npr1_40_list.txt'};
in_list = {'npr1_40_g_list.txt'}
using_second = 1;

% Second data set
%in_list = {'N2_40_g_list.txt', 'npr1_40_g_list.txt'};
%using_second = 1;

exp_ss_array = cell(length(in_list), num_statistics);
strain_vars = zeros(length(in_list),num_statistics-1);

% For each one of these list files provided i.e. for each group of videos..
for strain = 1:length(in_list)
    
    % --------------------- Finding the right files --------------------- %
    
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
        m
        
        % Manipulate the stored string to get the filename
        full = exp_file_locations{m};
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
        
        firstFrame = double(round(max(track_data.frame_number)*0.1));
        
        stationary_filter = frames < lastFrame & frames > firstFrame;
        x_data = x_data(stationary_filter);
        y_data = y_data(stationary_filter);
        frames = frames(stationary_filter);
        
        % Need to convert experimental data to simulation output structure
        % (TODO) ~ or add support to existing funtion
        in_data = {x_data, y_data, frames};
        
        % Then compute all chosen summary statistics, as with the simulated
        % data above.
        
        try
            % Compute SS1
            [SS1,SS2, SS3] = inf_clusterdistribution(in_data, 'experiment');
            exp_replicate_ss_array{m,2} = SS1;
            
            % Compute SS2
            exp_replicate_ss_array{m,3} = SS2;
            
            % Compute SS3
            exp_replicate_ss_array{m,4} = SS3;
            
        catch
            fprintf('failed SS1/SS2, but caught')
            exp_replicate_ss_array{m,2} = randi(10,1,40);
            exp_replicate_ss_array{m,3} = randi(10,1,40);
            exp_replicate_ss_array{m,4} = randi(10)
        end
        
        try
            % Continue for all chosen summary statistics
            % SS4: gr array, showing distribution of neighbours with distance
            exp_replicate_ss_array{m,5} = inf_pcf(in_data, 'experiment');
            
        catch
            fprintf('failed SS3, but caught')
            
            exp_replicate_ss_array{m,5} = randi(10,1,40);
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
















%*   Then, compute the appropriate distances between each of the
%    simple simulations and the experimental references

expsimple_dists = zeros(length(in_list),length(simple_sim_file_names), num_statistics);

for i = 1:length(in_list)
    for j = 1:length(simple_sim_file_names)
        
        if i ==1 & j ==1
            running_dists = zeros(1,size(exp_ss_array,2)-1,length(in_list));
        end
        
        exp_data = exp_ss_array(i,2:end);
        sim_data = simple_sim_ss_array(j,2:end);
        
        for ss = 1:length(sim_data)
            if length(exp_data{ss})>1
                exp_data{ss} = norm(exp_data{ss}-sim_data{ss});
                sim_data{ss} = 0;
            end
        end
        
        % Compute the distance between this simulation and the reference
        dists = cell2mat(exp_data)-cell2mat(sim_data);
        running_dists(:,:,i) = running_dists(:,:,i)+dists;
        expsimple_dists(i,j,:) = horzcat(norm(dists),dists);
    end
end

% Divide the distances from each ss by the mean for that strain
for strain = 1:length(in_list)
    for ss = 2:num_statistics
        strain_dist_means(strain, ss-1) = mean(expsimple_dists(strain,:,ss));
        expsimple_dists(strain,:,ss) = expsimple_dists(strain,:,ss)...
            ./mean(expsimple_dists(strain,:,ss));
    end
    
    for sim = 1:length(simple_sim_file_names)
        expsimple_dists(strain,sim,1) = sum(expsimple_dists(strain,sim,2:num_statistics));
    end
end















%*  Then, compute the appropriate distances between each of the
%   complex simulations and the experimental references

expcomplex_dists = zeros(length(in_list),length(complex_sim_file_names), num_statistics);

for i = 1:length(in_list)
    for j = 1:length(complex_sim_file_names)
        
        if i ==1 & j ==1
            running_dists = zeros(1,size(exp_ss_array,2)-1,length(in_list));
        end
        
        exp_data = exp_ss_array(i,2:end);
        sim_data = complex_sim_ss_array(j,2:end);
        
        for ss = 1:length(sim_data)
            if length(exp_data{ss})>1
                exp_data{ss} = norm(exp_data{ss}-sim_data{ss});
                sim_data{ss} = 0;
            end
        end
        
        % Compute the distance between this simulation and the reference
        dists = cell2mat(exp_data)-cell2mat(sim_data);
        running_dists(:,:,i) = running_dists(:,:,i)+dists;
        expcomplex_dists(i,j,:) = horzcat(norm(dists),dists);
    end
end

% Divide the distances from each ss by the mean for that strain
for strain = 1:length(in_list)
    for ss = 2:num_statistics
        strain_dist_means(strain, ss-1) = mean(expcomplex_dists(strain,:,ss));
        expcomplex_dists(strain,:,ss) = expcomplex_dists(strain,:,ss)...
            ./mean(expcomplex_dists(strain,:,ss));
    end
    
    for sim = 1:length(complex_sim_file_names)
        expcomplex_dists(strain,sim,1) = sum(expcomplex_dists(strain,sim,2:num_statistics));
    end
end













%* Then obtain the rough scaling map from simple to complex

mapping_parameters = {'revRateClusterEdge', 'vs'};
list_parameter_values_used = zeros(length(complex_sim_file_names),length(mapping_parameters),2);

for mapping_point = 1:length(complex_sim_file_names)
    
    load(simple_sim_file_names{mapping_point})
    for par = 1:length(mapping_parameters)
        list_parameter_values_used(mapping_point,par,1) = ...
            eval(strcat('param.', mapping_parameters{par}));
    end
    
    load(complex_sim_file_names{mapping_point})
    for par = 1:length(mapping_parameters)
        list_parameter_values_used(mapping_point,par,2) = ...
            eval(strcat('param.', mapping_parameters{par}));
    end
end

% Check to see if the parameters for the complex and simple simulations
% are identical, from the order inputted when files were read
if unique(list_parameter_values_used(:,:,1) == list_parameter_values_used(:,:,2)) ~= 1
    report_error = 'Simple and complex simulations have not been inputted in the same order'
else
    
    % Obtain the ratios between the simple and experimental distances
    rough_scales = expcomplex_dists(:,:,1) ./ expsimple_dists(:,:,1);
    
    % Sort the parameter pairs first by vs, then by revRate
    parameter_pairs_raw = list_parameter_values_used(:,:,1);
    [parameter_pairs_sorted, sorting_index] = sortrows(parameter_pairs_raw, [1,2]);
    
    % Order to scale factors in line with the parameters
    rough_scales_ordered = rough_scales(sorting_index);
    
    % Convert linear, ordered list to matrix
    rough_map = flipud(reshape(rough_scales_ordered, ...
        [length(unique(parameter_pairs_sorted(:,2))) length(unique(parameter_pairs_sorted(:,1)))]));
    
    % Check that the sorting and reshaping corretly orders the parameter
    % space
    parameter_pairs_revRate = parameter_pairs_raw(:,1);
    param_map_revRate = flipud(reshape(parameter_pairs_revRate(sorting_index), ...
        [length(unique(parameter_pairs_sorted(:,2))) length(unique(parameter_pairs_sorted(:,1)))]));
    
    parameter_pairs_vs = parameter_pairs_raw(:,2);
    param_map_vs = flipud(reshape(parameter_pairs_vs(sorting_index), ...
        [length(unique(parameter_pairs_sorted(:,2))) length(unique(parameter_pairs_sorted(:,1)))]));
    
    % How many intervals are wanted between cells in the rough map
    num_smoothing_intervals = 2;
    
    %Smooth the rough map with interpolation
    smooth_map = interp2(rough_map, num_smoothing_intervals);
    
    % Also smooth the corresponding parameter maps
    smooth_param_map_revRate = interp2(param_map_revRate, num_smoothing_intervals);
    smooth_param_map_vs = interp2(param_map_vs, num_smoothing_intervals);
    
    smooth_revRate_used = unique(smooth_param_map_revRate);
    smooth_vs_used = unique(smooth_param_map_vs);
end














%* Then read in the full cohort of simple simulations, and obtain distances
% to the same reference

if sim_ss_defined == 1
    full_sim_ss_array = sim_ss_array;
    full_sim_file_names = sim_file_names;
    
    full_sim_params = zeros(length(full_sim_file_names), length(mapping_parameters));
    
    % For each simulation file in the list, obtain the parameter values
    % used
    
    for sim = 1:length(full_sim_file_names)
        progress_reading_full_simulation_cohort = sim/length(full_sim_file_names)
        load(full_sim_file_names{sim});
        
        %Store the parameter values for later
        for par = 1:length(mapping_parameters)
            full_sim_params(sim,par) = ...
                eval(strcat('param.', mapping_parameters{par}));
        end
    end
else
    % Firstly, read in the list of simulation files
    full_sim_file_list = 'full_simulation_list.txt';
    
    % For speed now only !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    %full_sim_file_list = 'simple_simulation_list_mapping.txt';
    
    %sim_file_list = 'linus_sim_list.txt';
    full_sim_file_names = {};
    my_list_file = fopen(full_sim_file_list);
    
    % Read in the first line of the .txt file
    tline = fgetl(my_list_file);
    % For every line that still contains text...
    while ischar(tline)
        % Store this file location in the initialised array
        full_sim_file_names{end+1} = tline;
        % Read in the next line
        tline = fgetl(my_list_file);
    end
    % Close the list file, having read in all of the file names
    fclose(my_list_file)
    
    
    % Initialise an array for storing the summary statistic outputs
    num_statistics = 5;
    full_sim_ss_array = cell(length(full_sim_file_names), num_statistics);
    full_sim_params = zeros(length(full_sim_file_names), length(mapping_parameters));
    
    % For each simulation file in the list, compute the appropriate summary
    % statistics using supplied functions
    for sim = 1:length(full_sim_file_names)
        progress_reading_all_10000_simulations = sim/length(full_sim_file_names)
        load(full_sim_file_names{sim});
        
        %Store the parameter values for later
        for par = 1:length(mapping_parameters)
            full_sim_params(sim,par) = ...
                eval(strcat('param.', mapping_parameters{par}));
        end
        
        full_sim_ss_array{sim,1} = full_sim_file_names{sim};
        % Read in the appropriate data
        data = xyarray;
        % Compute SS1: Distribution of the number of worms per cluster
        [SS1,SS2,SS3] = inf_clusterdistribution(data, 'simulation');
        full_sim_ss_array{sim,2} = SS1;
        
        % SS2: Distribution of the branch heights from the same clustering
        full_sim_ss_array{sim,3} = SS2;
        
        % SS3 : Percentage of worms in cluster
        full_sim_ss_array{sim,4} = SS3;
        
        % SS4: gr array, showing distribution of neighbours with distance
        full_sim_ss_array{sim,5} = inf_pcf(data, 'simulation');
        
    end
end

expfullsimple_dists = zeros(length(in_list),length(full_sim_file_names), num_statistics);

for i = 1:length(in_list)
    for j = 1:length(full_sim_file_names)
        
        if i ==1 & j ==1
            running_dists = zeros(1,size(exp_ss_array,2)-1,length(in_list));
        end
        
        exp_data = exp_ss_array(i,2:end);
        sim_data = full_sim_ss_array(j,2:end);
        
        for ss = 1:length(sim_data)
            if length(exp_data{ss})>1
                exp_data{ss} = norm(exp_data{ss}-sim_data{ss});
                sim_data{ss} = 0;
            end
        end
        
        % Compute the distance between this simulation and the reference
        dists = cell2mat(exp_data)-cell2mat(sim_data);
        running_dists(:,:,i) = running_dists(:,:,i)+dists;
        expfullsimple_dists(i,j,:) = horzcat(norm(dists),dists);
    end
end

% Divide the distances from each ss by the mean for that strain
for strain = 1:length(in_list)
    for ss = 2:num_statistics
        strain_dist_means(strain, ss-1) = mean(expfullsimple_dists(strain,:,ss));
        expfullsimple_dists(strain,:,ss) = expfullsimple_dists(strain,:,ss)...
            ./mean(expfullsimple_dists(strain,:,ss));
    end
    
    for sim = 1:length(full_sim_file_names)
        expfullsimple_dists(strain,sim,1) = sum(expfullsimple_dists(strain,sim,2:num_statistics));
    end
end











%* Scale the 10,000 simulations by the new fine map
projected_complex_dists = expfullsimple_dists(:,:,1)

for i = 1:length(projected_complex_dists(:,:,1))
    
    %Recover the value of revRate and vs
    revRate_value = full_sim_params(i,1);
    vs_value = full_sim_params(i,2);
    
    [M,revRate_map_index] = min(abs(smooth_revRate_used(:)-revRate_value))
    [M,vs_map_index] = min(abs(smooth_vs_used(:)-vs_value))
    
    dist_scaler = smooth_map(vs_map_index, revRate_map_index);
    
    projected_complex_dists(1,i,1) = projected_complex_dists(1,i,1) * dist_scaler;
    
end









%% Plotting posteriors
%  Produce new posterior by considering the new distances in the
%  more complex simulation space

% We want to accept the best n% of simulations
% Set the cutoffs for taking the top n% of simulations
n_cuts = [0.01];
test_params = {'vs', 'revRateClusterEdge'};

chosen_params = zeros(floor(prod(size(projected_complex_dists))*max(n_cuts)/num_statistics),...
    length(test_params), length(n_cuts));

% For each of these cutoffs, produce distributions of the parameters
for cutoff = 1:length(n_cuts)
    n = n_cuts(cutoff);
    top_n = floor(length(projected_complex_dists))*n;
    
    lin = reshape(projected_complex_dists(:,:,1).' ,1,numel(projected_complex_dists(:,:,1)));
    A = sort(lin);
    %A = sort(expsim_dists(:,:,1));
    B = (projected_complex_dists(:,:,1)<=A(top_n)).*projected_complex_dists(:,:,1);
    
    best_sims = floor(find(B)/length(in_list));
    list_best = {};
    
    for sim = 1:length(best_sims)
        list_best{end+1} = full_sim_file_names(best_sims(sim));
    end
    
    % Use the mat file to find the parameters for plotting, no need to reload
    % each of the simulations
    
    for par = 1:length(test_params)
        
        for i = 1:length(list_best)
            load(list_best{i}{1});
            chosen_params(i,par,cutoff) = eval(strcat('param.', test_params{par}));
        end
        
    end
    
end

% Comparing joint distributions of parameters
figure;
for cutoff = 1:length(n_cuts)
    to_plot = chosen_params(:,:,cutoff);
    
    %Eliminate redundantrows, where all parameter values are zero
    to_plot = to_plot(any(to_plot~=0,2),:);
    
    subplot(1,length(n_cuts),cutoff)
    [S,AX,BigAx,H,HAx] = plotmatrix(to_plot);
    
    title(['Map adjusted prior, using top ' num2str(n_cuts(cutoff)*100) '% of simulations'])
    
    for i = 1:length(test_params)
        ylabel(AX(i,1),test_params(i))
        xlabel(AX(length(test_params),i),test_params(i))
    end
    
end





% And repeat for the original simple simulations, for comparison

% We want to accept the best n% of simulations
% Set the cutoffs for taking the top n% of simulations
n_cuts = [0.01];



chosen_params = zeros(floor(prod(size(expfullsimple_dists))*max(n_cuts)/num_statistics),...
    length(test_params), length(n_cuts));


% For each of these cutoffs, produce distributions of the parameters
for cutoff = 1:length(n_cuts)
    n = n_cuts(cutoff)
    top_n = floor(length(projected_complex_dists))*n;
    
    lin = reshape( expfullsimple_dists(:,:,1).' ,1,numel(expfullsimple_dists(:,:,1)));
    A = sort(lin);
    %A = sort(expsim_dists(:,:,1));
    B = (expfullsimple_dists(:,:,1)<=A(top_n)).*expfullsimple_dists(:,:,1);
    
    best_sims = floor(find(B)/length(in_list));
    list_best = {};
    
    for sim = 1:length(best_sims)
        list_best{end+1} = full_sim_file_names(best_sims(sim));
    end
    
    % Use the mat file to find the parameters for plotting, no need to reload
    % each of the simulations
    
    for par = 1:length(test_params)
        
        for i = 1:length(list_best)
            load(list_best{i}{1});
            chosen_params(i,par,cutoff) = eval(strcat('param.', test_params{par}));
        end
        
    end
    
end

% Comparing joint distributions of parameters
figure;
for cutoff = 1:length(n_cuts)
    to_plot = chosen_params(:,:,cutoff);
    
    %Eliminate redundantrows, where all parameter values are zero
    to_plot = to_plot(any(to_plot~=0,2),:);
    
    subplot(1,length(n_cuts),cutoff)
    [S,AX,BigAx,H,HAx] = plotmatrix(to_plot);
    
    title(['Unadjusted prior, using top ' num2str(n_cuts(cutoff)*100) '% of simulations'])
    
    for i = 1:length(test_params)
        ylabel(AX(i,1),test_params(i))
        xlabel(AX(length(test_params),i),test_params(i))
    end
    
end