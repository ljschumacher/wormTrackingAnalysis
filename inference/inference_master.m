% Setting parameters
% State how many summary statistics are to be computed
num_statistics = 1

addpath('component_functions')

%% Analyse simulation data
sim_ss_array = f_analyse_sims('datalists/woidM18_2750ish_samples.txt', 0, num_statistics);

%% Analyse experimental data
[exp_ss_array, exp_strain_list] = f_analyse_exps(...
    {'npr1'},2,num_statistics);

%% Obtain distances between each of the experiments and simulations
expsim_dists = f_exp2sim_dist(...
    exp_ss_array, sim_ss_array, exp_strain_list);

%% Perform parameter inference
[chosen_params, chosen_samples] = f_infer_params(...
    expsim_dists, {'dkdN', 'revRateClusterEdge'}, [0.01],...
    '../../../sworm-model/woidModel/paramSamples_nSim10000_nParam2.mat');
