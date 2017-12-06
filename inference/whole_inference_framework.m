%% C.elegans parameter inference framework
% This is a framework for approximating distributions of model parameters
% by comparing the distances between experimental data and simulations of
% worm behaviour. 

% Critically, a list of the experimental recordings to be used as
% references must be provided and assigned to the variable 'in_list',
% alongside a list of the simulations - assigned to 'sim_file_list'.

% The files named in these lists must be present in the current working
% directory.

% This framework is divided into a number of sections for ease of use:
% --- 1. Analysing simulation data
% --- 2. Building the experimental reference
% --- 3. Computing distances between experiments and simulations
% --- 4. Accepting the closest n% of simulations
% --- 5. Producing distributions of inferred parameters

% OPTIONAL
% --- 6. Reporting further details


%% ---------------------- Analysing simulation data ---------------------- % 
close all; clear; clc;

% Firstly, read in the list of simulation files
sim_file_list = 'list_extended_simulations.txt';
sim_file_names = {};
my_list_file = fopen(sim_file_list);

% Read in the first line of the .txt file
tline = fgetl(my_list_file);
% For every line that still contains text...
while ischar(tline)
    % Store this file location in the initialised array
    sim_file_names{end+1} = tline;
    % Read in the next line
    tline = fgetl(my_list_file);
end
% Close the list file, having read in all of the file names
fclose(my_list_file)

% State how many summary statistics are to be computed
num_statistics = 4;


% Construct array for storing the outputs of the summary statistics for
% each simulation. Firstly, have to account for the extra column for the 
% strain/exp name
num_statistics = num_statistics+1;
sim_ss_array = cell(length(sim_file_names), num_statistics);

% For each simulation file in the list, compute the appropriate summary
% statistics using supplied functions
for sim = 1:length(sim_file_names)
    progress = sim/length(sim_file_names)
    load(sim_file_names{sim});
    sim_ss_array{sim,1} = sim_file_names{sim};
    % Read in the appropriate data
    data = xyarray;
    % Compute SS1: Distribution of the number of worms per cluster
    [SS1,SS2,SS3] = inf_clusterdistribution(data, 'simulation');
    sim_ss_array{sim,2} = SS1;
    
    % SS2: Distribution of the branch heights from the same clustering
    sim_ss_array{sim,3} = SS2;
    
    % SS3 : Percentage of worms in cluster
    sim_ss_array{sim,4} = SS3;
    
    % SS4: gr array, showing distribution of neighbours with distance
    sim_ss_array{sim,5} = inf_gr(data, 'simulation');
    
    % Continue for all chosen summary statistics
    
end

%% ----------------- Building the experimental reference -----------------%
% Then consider the experimental data, which will become the references for
% computing distances to simulations.

% Supply the list files, which contain the directory locations of the
% videos for the titular strains. Should specify whether data is from the
% first or second set of experiments with 'using_second = 0' and
% 'using_second = 1' respectively. This changes the minimum brightness used
% to filter out larvae.

% First data set
in_list = {'N2_40_list.txt', 'HA_40_list.txt', 'npr1_40_list.txt'};
in_list = {'N2_40_g_list.txt'}
using_second = 1;

% Construct arrays for storing the summary statistic outputs
exp_ss_array = cell(length(in_list), num_statistics);
strain_vars = zeros(length(in_list),num_statistics-1);

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
        
        try
            % Compute SS1
            [SS1,SS2, SS3] = inf_clusterdistribution(in_data, 'experiment');
            exp_replicate_ss_array{m,2} = SS1;
            
            % Compute SS2
            exp_replicate_ss_array{m,3} = SS2;
            
            % Compute SS3
            exp_replicate_ss_array{m,4} = SS3;
            
        catch
            fprintf('failed SS1/SS2/SS3, but caught')
            exp_replicate_ss_array{m,2} = randi(10,1,40);
            exp_replicate_ss_array{m,3} = randi(10,1,40);
            exp_replicate_ss_array{m,4} = randi(10)
        end
        
        try
            % Continue for all chosen summary statistics
            % SS4: gr array, showing distribution of neighbours with distance
            exp_replicate_ss_array{m,5} = inf_gr(in_data, 'experiment');
            
        catch
            fprintf('failed SS4, but caught')
            
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


%% ------- Computing distances between experiments and simulations ------ %
% Then, compute the appropriate distances between each of the
% simulations and the experimental references

% Constuct a n-by-m-by-l+1 matrix for containing the distances between each 
% experiment (n) and each simulation (m) for each of the summary statistics
% computed (l)

expsim_dists = zeros(length(in_list),length(sim_file_names), num_statistics);

for i = 1:length(in_list)
    for j = 1:length(sim_file_names)
        
        if i ==1 & j ==1
            running_dists = zeros(1,size(exp_ss_array,2)-1,length(in_list));
        end
        
        exp_data = exp_ss_array(i,2:end);
        sim_data = sim_ss_array(j,2:end);
        
        for ss = 1:length(sim_data)
            if length(exp_data{ss})>1
                exp_data{ss} = norm(exp_data{ss}-sim_data{ss});
                sim_data{ss} = 0;
            end
        end
        
        % Compute the distance between this simulation and the reference
        dists = cell2mat(exp_data)-cell2mat(sim_data);
        running_dists(:,:,i) = running_dists(:,:,i)+dists;
        expsim_dists(i,j,:) = horzcat(norm(dists),dists);
    end
end

% Normalization: divide the distances from each ss by the mean for that
% strain. This ensures that each summary statistic contributes is weighted
% similarly in the distance function. Else, summary statistics at a larger
% scale i.e. % of worms in cluster, would dwarf others when the mean
% euclidean distance is calculated.

for strain = 1:length(in_list)
    for ss = 2:num_statistics
        strain_dist_means(strain, ss-1) = mean(expsim_dists(strain,:,ss));
        expsim_dists(strain,:,ss) = expsim_dists(strain,:,ss)...
            ./mean(expsim_dists(strain,:,ss));  
    end
    
    for sim = 1:length(sim_file_names)
        expsim_dists(strain,sim,1) = sum(expsim_dists(strain,sim,2:num_statistics));
    end
end

% % OR: divide the distances from each ss by the mean across all strains
%     for ss = 2:num_statistics
%         expsim_dists(:,:,ss) = expsim_dists(:,:,ss)...
%             ./mean(expsim_dists(:,:,ss));  
%     end
%     
% for strain = 1:length(in_list)    
%     for sim = 1:length(sim_file_names)
%         expsim_dists(strain,sim,1) = sum(expsim_dists(strain,sim,2:num_statistics));
%     end
% end


% Consider the distance composition to check whether this normalization of
% the summary statistic weights has been successful.
sim_strain_vars = zeros(length(in_list),num_statistics-1);
strain_dist_means = zeros(length(in_list),num_statistics-1);

for strain = 1:length(in_list)
    for ss = 2:num_statistics
        sim_strain_vars(strain, ss-1) = var(expsim_dists(strain,:,ss));
    end
end

% To ensure pie is full
if sum(sim_strain_vars)<=1
    sim_strain_vars(:) = sim_strain_vars(:).*(1/sum(sim_strain_vars));
end

pie_data = zeros(strain,num_statistics-1);
figure;
for strain = 1:length(in_list)
    for ss = 1:num_statistics-1
        pie_data(strain,ss) = sum(expsim_dists(strain,:,ss+1));
    end
        
    subplot(2,length(in_list)+1,strain)
    pie(pie_data(strain,:))
    title(in_list(strain), 'interpreter','none')
    
    subplot(2,length(in_list)+1,strain+length(in_list)+1)
    pie(sim_strain_vars(strain,:))
    title('var in distances')
end
subplot(2,length(in_list)+1,strain+1)
pie(sum(pie_data))
title('Average across strains')



%% -------------- Accepting the closest n% of simulations --------------- % 

% Set the cutoffs for taking the top n% of simulations e.g to select the
% closest 1% of simulations, use 'n_cuts = [0.01]'. To see the effect that
% using different cutoffs has on the parameter distributions inferred,
% separate values with a comma: 'n_cuts = [0.10,0.05,0.01].

n_cuts = [0.01];
test_params = {'vs', 'revRateClusterEdge', 'Rir', 'Ris'};

%Create array for storing the parameters of the top n% of simulations
chosen_params = zeros(floor(prod(size(expsim_dists))*max(n_cuts)/num_statistics),...
    length(test_params), length(n_cuts));

% For each of the % cutoffs specified in n_cuts, produce distributions of 
% the parameters
for cutoff = 1:length(n_cuts)
    n = n_cuts(cutoff)
    top_n = floor(prod(size(expsim_dists))*n/num_statistics);
    
    lin = reshape( expsim_dists(:,:,1).' ,1,numel(expsim_dists(:,:,1)));
    A = sort(lin);
    %A = sort(expsim_dists(:,:,1));
    B = (expsim_dists(:,:,1)<=A(top_n)).*expsim_dists(:,:,1);
    
    best_sims = floor(find(B)/length(in_list));
    list_best = {};
    
    for sim = 1:length(best_sims)
        list_best{end+1} = sim_file_names(best_sims(sim)+1);
    end
    
    % Could use the mat file to find the parameters for plotting, 
    % no need to reload each of the simulations. Remains fast without.
    
    for par = 1:length(test_params)
        for i = 1:length(list_best)
            load(list_best{i}{1});
            chosen_params(i,par,cutoff) = eval(strcat('param.', test_params{par}));
        end  
    end    
end


%% -------- Producing joint distributions of inferred parameters -------- %

figure;
for cutoff = 1:length(n_cuts)
    to_plot = chosen_params(:,:,cutoff);
    
    % Eliminate redundantrows, where all parameter values are zero
    % Useful when there are multiple, increasingly tight cutoffs
    to_plot = to_plot(any(to_plot~=0,2),:);
    
    subplot(1,length(n_cuts),cutoff)
    [S,AX,BigAx,H,HAx] = plotmatrix(to_plot);
    
    title(['Top ' num2str(n_cuts(cutoff)*100) '% of simulations'])
    
    for i = 1:length(test_params)
        ylabel(AX(i,1),test_params(i))
        xlabel(AX(length(test_params),i),test_params(i))
    end
    
end

%% --------------- Reporting further details [optional] ----------------- %
% From the top n% of simulations, it is possible to report information
% regarding the single closest simulation and its parameters. 

top_params = chosen_params(:,:,size(chosen_params,3)); 

%Eliminate redundantrows, where all parameter values are zero
top_params = top_params(any(top_params~=0,2),:);
[chosen_sim_dists, sort_index] = sort(B(B~=0));
chosen_sim_dists = chosen_sim_dists ./ min(chosen_sim_dists);
chosen_sim_dists = 1 ./ (chosen_sim_dists);

for p = 1:size(top_params,2)
to_sort = top_params(:,p);
top_params(:,p) = to_sort(sort_index);
end

% Report the closest single simulation match
list_best = list_best(sort_index)
single_best_match = list_best{1};
single_best_match = single_best_match{1};
single_best_params = top_params(1,:)

% Also, produce a weighted average of each parameter across the top n% of
% simulations. These are weighted so that the parameters of simulations
% twice as far from the experiments as the single best simulation are
% weighted by half in the computation of the weighted parameter mean.
pure_means = mean(top_params);
weighted_params = top_params;

for i = 1:length(chosen_sim_dists)
    weighted_params(i,:) = top_params(i,:).*chosen_sim_dists(i);
end

weighted_means = sum(weighted_params)/sum(chosen_sim_dists)