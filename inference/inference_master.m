% Setting parameters
% State how many summary statistics are to be computed

% issues/to-do:
% - calculation of summary statistics could be sped up by calculating all
% stats within the loop over frames, rather than looping over frames for
% each stat
num_statistics = 4
addpath('component_functions');

%% Analyse simulation data
sim_file_list = 'datalists/woidM18_20k_samples_N2like.txt';
filepath = '../../../sworm-model/woidModel/results/paramSampleResults/woidlinos/woidM18paramD2/';

[sim_ss_array, sim_file_names, param_return] = f_analyse_sims(sim_file_list,...
    filepath, {'revRateClusterEdge','dkdN_dwell'}, num_statistics);

%% Analyse experimental data
[exp_ss_array, exp_strain_list] = f_analyse_exps(...
    {'npr1','N2'},2,num_statistics);

%% Obtain distances between each of the experiments and simulations
expsim_dists = f_exp2sim_dist(exp_ss_array, sim_ss_array,[0.01, 0.01, 0.01, 10]);

%% optimise weightings of summary statistics
expsim_dists = f_optim_posterior(exp_ss_array(1,:), sim_ss_array,...
    0.025, '../../../sworm-model/woidModel/paramSamples_nSim20000_nParam2.mat');
%% Perform parameter inference
load('../../../sworm-model/woidModel/paramSamples_nSim20000_nParam2.mat')
[chosen_params, chosen_samples] = f_infer_params(...
    expsim_dists, exp_strain_list,[0.0125],paramSamples, true);

%% Plot summary statistics of experiments and best samples
for statCtr = 1:num_statistics
    figure
    for strainCtr = 1:length(exp_strain_list)
    plot(exp_ss_array{strainCtr,statCtr+1},'LineWidth',2)
    hold on
    end
        for strainCtr = 1:length(exp_strain_list)
    for ii=1:10
        plot(sim_ss_array{chosen_samples(strainCtr,ii),statCtr+1})
    end
        end
    title(['summary statistic ' num2str(statCtr)])
    legend([exp_strain_list{1} ' mean'],[exp_strain_list{2} ' mean'],'best simulations')
end

% %% plot surface of dissimilatirity
% load ../../../sworm-model/woidModel/paramSamples_nSim10000_nParam2.mat
% F = scatteredInterpolant(paramSamples.dkdN,paramSamples.revRateClusterEdge,expsim_dists(:,1));