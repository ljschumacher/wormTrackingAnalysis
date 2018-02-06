% Declare the two parameters for which to build to 2D map, as a cell array
mapping_params = {'revRateClusterEdge','dkdN_dwell'};

num_statistics = 3
addpath('component_functions');
sum_stat_weights = [0.27 0.71 0.02];
%% Analyse a small set of simple simulations first
filepath = '../../../sworm-model/woidModel/results/woidlinos/mapping/';
[simple_sim_ss_array, simple_sim_file_names, simple_param_return] = f_analyse_sims(...
    'datalists/wlM18_mapping_sims1.txt',filepath,mapping_params,num_statistics);

%% Analyse the corresponding complex simulations
filepath = '../../../sworm-model/woidModel/results/woids/mapping/';
[complex_sim_ss_array, complex_sim_file_names, complex_param_return] = f_analyse_sims(...
    'datalists/wM36_mapping_sims1.txt',filepath,mapping_params,num_statistics,'complexsim');

%% Then build/load the experimental reference with which to compare 
[exp_ss_array, exp_strain_list] = f_analyse_exps(...
    {'npr1','N2'},2,num_statistics);

%% Compute distances between the simple simulations and the experiments
expsimple_dists = f_exp2sim_dist(exp_ss_array,simple_sim_ss_array,sum_stat_weights);

%% Then compute the distances between the complex simulations and the exps
expcomplex_dists = f_exp2sim_dist(exp_ss_array,complex_sim_ss_array,sum_stat_weights);

%% Create a mapping
[map, smooth_p1, smooth_p2] = f_create_map(simple_param_return, complex_param_return, ...
    expsimple_dists, expcomplex_dists, mapping_params);

%% Analyse/load the full cohort of simple simulations
sim_file_list = 'datalists/woidM18_10k_samples_npr1like_reducedPrior.txt';
filepath = '../../../sworm-model/woidModel/results/paramSampleResults/woidlinos/woidM18paramD2/';

[sim_ss_array, sim_file_names, param_return] = f_analyse_sims(sim_file_list,...
    filepath, mapping_params, num_statistics);

%% Compute the initial distance from these simple simulations to the exps
expsim_dists = f_exp2sim_dist(exp_ss_array, sim_ss_array,sum_stat_weights);
    
%% Scale the distances for the full cohort of simple simulations by the map
projected_dists = f_apply_map( ...
    expsim_dists, param_return, map, smooth_p1, smooth_p2);

%% Perform parameter inference
load('../../../sworm-model/woidModel/paramSamples_nSim10079_nParam2_reducedPrior.mat')
[chosen_params, chosen_samples] = f_infer_params(...
    projected_dists, exp_strain_list,[0.125],paramSamples, true);
