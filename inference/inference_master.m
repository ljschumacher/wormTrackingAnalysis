% Setting parameters
% State how many summary statistics are to be computed
global num_statistics
num_statistics = 1
num_statistics = num_statistics+1;

addpath('component_functions')

%% Analyse simulation data
sim_ss_array = f_analyse_sims('datalists/woidM18_999samples.txt', 0);

%% Analyse experimental data
[exp_ss_array, exp_strain_list] = f_analyse_exps(...
    {'npr1'},2);

%% Obtain distances between each of the experiments and simulations
expsim_dists = f_exp2sim_dist(...
    exp_ss_array, sim_ss_array, exp_strain_list);

%% Perform parameter inference
chosen_params = f_infer_params(...
    expsim_dists, {'dkdN_dwell', 'revRateClusterEdge'}, [0.01]);
