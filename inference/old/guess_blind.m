%% Framework for testing whether parameters of blind simulations can be recovered
%  experimental data

close all; clear; clc;

% Firstly, read in the list of simulation files
sim_file_list = 'full_simulation_list.txt';
%sim_file_list = 'linus_sim_list.txt';
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


% Initialise an array for storing the summary statistic outputs
num_statistics = 5;
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
    sim_ss_array{sim,4} = SS3
    
    % SS4: gr array, showing distribution of neighbours with distance
    sim_ss_array{sim,5} = inf_gr(data, 'simulation');
    
    % Continue for all chosen summary statistics
    
end

%% Building the blind references
% Then consider the blind simulation data, which will become the reference
% for computing distances to simulations.

% Supply the list files, which contain the directory locations of the
% videos for the titular strains


% Firstly, read in the list of blind simulation files
blind_file_list = 'blind_trials.txt';
blind_file_list = 'full_blind_list.txt';

blind_sim_file_names = {};
my_list_file = fopen(blind_file_list);

% Read in the first line of the .txt file
tline = fgetl(my_list_file);
% For every line that still contains text...
while ischar(tline)
    % Store this file location in the initialised array
    blind_sim_file_names{end+1} = tline;
    % Read in the next line
    tline = fgetl(my_list_file);
end
% Close the list file, having read in all of the file names
fclose(my_list_file)


% Initialise an array for storing the summary statistic outputs
blind_sim_ss_array = cell(length(blind_sim_file_names), num_statistics);

% For each simulation file in the list, compute the appropriate summary
% statistics using supplied functions
for sim = 1:length(blind_sim_file_names)
    progress = sim/length(blind_sim_file_names)
    load(blind_sim_file_names{sim});
    blind_sim_ss_array{sim,1} = blind_sim_file_names{sim};
    % Read in the appropriate data
    data = xyarray;
    % Compute SS1: Distribution of the number of worms per cluster
    [SS1,SS2,SS3] = inf_clusterdistribution(data, 'simulation');
    blind_sim_ss_array{sim,2} = SS1;
    
    % SS2: Distribution of the branch heights from the same clustering
    blind_sim_ss_array{sim,3} = SS2;
    
    % SS3 : Percentage of worms in cluster
    blind_sim_ss_array{sim,4} = SS3;
    
    % SS4: gr array, showing distribution of neighbours with distance
    blind_sim_ss_array{sim,5} = inf_gr(data, 'simulation');
    
    % Continue for all chosen summary statistics
end


% COMPUTING DISTANCES
% Then, compute the appropriate distances between each of the simulations 
% and the blind references
exp_ss_array = blind_sim_ss_array;
expsim_dists = zeros(length(blind_sim_file_names),length(sim_file_names), num_statistics);

for i = 1:length(blind_sim_file_names)
    
    for j = 1:length(sim_ss_array) 
        
        
        if i ==1 & j ==1
            running_dists = zeros(1,size(exp_ss_array,2)-1,length(blind_sim_file_names));
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

% NORMALISATION
% Divide the distances from each ss by the mean for that strain
for strain = 1:length(blind_sim_file_names)
    for ss = 2:num_statistics
        strain_dist_means(strain, ss-1) = mean(expsim_dists(strain,:,ss));
        expsim_dists(strain,:,ss) = expsim_dists(strain,:,ss)...
            ./mean(expsim_dists(strain,:,ss));
    end
    
    for sim = 1:length(sim_ss_array) 
        expsim_dists(strain,sim,1) = sum(expsim_dists(strain,sim,2:num_statistics));
    end
end


% We want to accept the best n% of simulations
test_params = {'revRateClusterEdge','vs', 'Ris', 'Rir'};
hidden_params = load('blindSamples_nSim100_nParam4.mat');

% For storing the raw predicted and real parameter values for the blind
% simulations
parameter_store = zeros(length(blind_sim_file_names), length(test_params),2);

% For storing the log ratio of real to predicted parameters, along side the
% confidence interval the real paramter falls into, based on the posterioir
% distribution of the prediction
predictions_to_real = zeros(length(blind_sim_file_names), length(test_params),2);

figure;

for blind_trial = 1:length(blind_sim_file_names)
    blind_trial
    % Set the cutoffs for taking the top n% of simulations
    n_cuts = [0.01];
    
    chosen_params = zeros(floor(prod(size(expsim_dists))*max(n_cuts)/num_statistics),...
        length(test_params), length(n_cuts));
    
    % For each of these cutoffs, produce distributions of the parameters
    for cutoff = 1:length(n_cuts)
        n = n_cuts(cutoff);
        top_n = size(expsim_dists,2)*n;
        
        lin = reshape(expsim_dists(blind_trial,:,1).' ,1,numel(expsim_dists(blind_trial,:,1)));
        A = sort(lin);
        %A = sort(expsim_dists(:,:,1));
        B = (expsim_dists(blind_trial,:,1)<=A(top_n)).*expsim_dists(blind_trial,:,1);
        
        %best_sims = floor(find(B)/length(blind_sim_file_names));
        best_sims = find(B);
        list_best = {};
        
        for sim = 1:length(best_sims)
            list_best{end+1} = sim_file_names(best_sims(sim));
        end
        
        % Use the mat file to find the parameters for plotting, no need to reload
        % each of the simulations
        
        for par = 1:length(test_params)
            
            for q = 1:length(list_best)
                load(list_best{q}{1});
                chosen_params(q,par,cutoff) = eval(strcat('param.', test_params{par}));
            end
            
        end
        
    end
 
    
    % Producing an estimate for the reference/blind-sim parameters
    top_params = chosen_params(:,:,size(chosen_params,3));
    
    %Eliminate redundantrows, where all parameter values are zero
    top_params = top_params(any(top_params~=0,2),:);
    
    [blind_sim_dists, sort_index] = sort(B(B~=0));
    blind_sim_dists = blind_sim_dists ./ min(blind_sim_dists);
    blind_sim_dists = 1 ./ (blind_sim_dists);
    
    for p = 1:size(top_params,2)
        to_sort = top_params(:,p);
        top_params(:,p) = to_sort(sort_index);
    end
    
    % Report the closest single simulation match
    list_best = list_best(sort_index);
    single_best_match = list_best{1};
    single_best_match = single_best_match{1};
    single_best_params = top_params(1,:);
    
    pure_means = mean(top_params);
    
    weighted_params = top_params;
    
    for i = 1:length(blind_sim_dists)
        weighted_params(i,:) = top_params(i,:).*blind_sim_dists(i);
    end
    
    weighted_means = sum(weighted_params)/sum(blind_sim_dists);
    
    % For each parameter being fitted, get the log ratio and CI
    
    for model_parameter = 1:length(test_params);
        predicted = weighted_means(model_parameter);
        real_value = hidden_params.paramSamples{blind_trial,model_parameter};
        difference = abs(predicted-real_value);
        
        % Store the raw parameter values
        parameter_store(blind_trial, model_parameter,1) = predicted;
        parameter_store(blind_trial, model_parameter,2) = real_value;
        
        
        %Compare the real value to the predicted value with a log ratio
        predictions_to_real(blind_trial, model_parameter,1) = ...
            abs(log(predicted/real_value));
       % Also report the associated confidence score
       to_search = top_params(:,model_parameter);
       to_search = to_search(to_search > (predicted-difference) & to_search < (predicted+difference));

        
       predictions_to_real(blind_trial, model_parameter,2) = ...
           1 - ((length(to_search))/length(top_params));
    end
end

figure;
for par = 1:length(test_params)
    subplot(2,2,par)
    scatter(parameter_store(:,par,1), parameter_store(:,par,2))
    hold on
    
    refline(1,0)
    xlabel('Predicted value')
    ylabel('Real value')
    title(test_params{par})
end
