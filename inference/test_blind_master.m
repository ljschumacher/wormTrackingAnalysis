%% Framework for testing whether parameters of blind simulations can be recovered
%  experimental data

num_statistics = 4

addpath('component_functions');

%% Analyse simulation data - or load precomputed summary stats
sim_file_list = 'datalists/woidM18_10k_samples.txt';
filepath = '../../../sworm-model/woidModel/results/paramSampleResults/woidlinos/woidM18paramD2/';

[sim_ss_array, sim_file_names, param_return] = f_analyse_sims(sim_file_list,...
    filepath, {'revRateClusterEdge','dkdN_dwell'}, num_statistics);

%% Building the blind references
% Then consider the blind simulation data, which will become the reference
% for computing distances to simulations.
blind_file_list = 'datalists/woidM18_100_blind_samples.txt';
filepath = '../../../sworm-model/woidModel/results/paramSampleResults/woidlinos/woidM18paramD2/blind/';

[blind_ss_array, blind_file_names] = f_analyse_sims(blind_file_list,filepath,0,num_statistics);

%% COMPUTING DISTANCES
% Then, compute the appropriate distances between each of the simulations
% and the blind references
blindsim_dists = f_exp2sim_dist(...
    blind_ss_array, sim_ss_array, blind_file_names);

%% infer parameters
test_params = {'revRateClusterEdge','dkdN'};
param_file = '../../../sworm-model/woidModel/paramSamples_nSim20000_nParam2.mat';
[chosen_params, chosen_samples] = f_infer_params(blindsim_dists,...
    blind_file_names,test_params, 0.01, param_file, false);

%% evaluate quality of inference
hidden_params = load('../../../sworm-model/woidModel/blindSamples_nSim100_nParam2.mat',...
    'paramSamples');

predictions_to_real = f_evaluate_inference(chosen_samples,...
    chosen_params,blindsim_dists,blind_file_names,hidden_params,test_params,true);

