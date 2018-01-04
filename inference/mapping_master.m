% Declare the two parameters for which to build to 2D map, as a cell array
mapping_params = {'revRateClusterEdge', 'vs'}

% State how many summary statistics will be computed (then increment by 1)
num_statistics = 4

% Analyse a small set of simple simulations first
[simple_sim_ss_array, simple_sim_file_names] = f_analyse_sims(...
    'simple_simulation_list_mapping.txt')

% Analyse the corresponding complex simulations
[complex_sim_ss_array, complex_sim_file_names] = f_analyse_sims(...
    'complex_simulation_list_mapping.txt')

% Then build the experimental reference with which to compare 
% the simple and complex simulations
[exp_ss_array, exp_strain_list] = f_analyse_exps(...
    {'N2_40_g_list.txt'},1)

% Compute distances between the simple simulations and the experiments
expsimple_dists = f_exp2sim_dist(...
    exp_ss_array, simple_sim_ss_array, exp_strain_list)

% Then compute the distances between the complex simulations and the exps
expcomplex_dists = f_exp2sim_dist(...
    exp_ss_array, complex_sim_ss_array, exp_strain_list)

% Create a mapping
[map, smooth_p1, smooth_p2] = f_create_map(simple_sim_file_names, complex_sim_file_names, ...
    expsimple_dists, expcomplex_dists, mapping_params)

% Analyse the full cohort of simple simulations
[full_simple_sim_ss_array, full_simple_sim_file_names, extracted_params]...
    = f_analyse_sims('simple_simulation_list_mapping.txt', mapping_params)

% Compute the initial distance from these simple simulations to the exps
full_expsimple_dists = f_exp2sim_dist(...
    exp_ss_array, full_simple_sim_ss_array, exp_strain_list)

    
% Scale the distances for the full cohort of simple simulations by the map
projected_dists = apply_dist_map( ...
    full_expsimple_dists, full_sim_params, smooth_map, smooth_p1, smooth_p2)


% Perform parameter inference
chosen_params = f_infer_params(...
    projected_dists, {'vs', 'revRateClusterEdge', 'Rir', 'Ris'}, [0.01])
