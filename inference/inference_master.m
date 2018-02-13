function [weights_optim, min_obj] = inference_master(model,strain,accept_ratio)
% Inference framework

% issues/to-do:
% - calculation of summary statistics could be sped up by calculating all
% stats within the loop over frames, rather than looping over frames for each stat
addpath('component_functions');

% % set overall parameters
% model = 'log-rods'; % 'rods' or 'worms'
% strain = 'npr1';
% accept_ratio = 0.02;
param_names = {'revRateClusterEdge','dkdN_dwell'};

switch model
    case 'rods'
        num_statistics = 3;
        load('../../../sworm-model/woidModel/paramSamples_nSim20000_nParam2.mat')
        sumstat_filename = ['sumstats_20ksamples_wlM18_' strain 'like.mat'];
        sim_file_list = ['datalists/woidM18_20k_samples_' strain 'like.txt'];
        filepath = '../../../sworm-model/woidModel/results/paramSampleResults/woidlinos/woidM18paramD2/';
        scaleflag = 'linear';
     case 'log-rods'
        num_statistics = 3;
        load('../../../sworm-model/woidModel/paramSamples_log_nSim30000_nParam2.mat')
        sumstat_filename = ['sumstats_30klogsamples_wlM18_' strain 'like.mat'];
        sim_file_list = ['datalists/woidM18_30k_logsamples_' strain 'like.txt'];
        filepath = '../../../sworm-model/woidModel/results/paramSampleResults/paramSamplesLog/woidlinos/';
        scaleflag = 'log';
    case 'worms'
        num_statistics = 3; % 4th stat, polar order, did not seem to work well
        load('../../../sworm-model/woidModel/paramSamples_nSim10000_nParam2.mat')
        sumstat_filename = ['sumstats_10ksamples_wM36_' strain 'like.mat'];
        sim_file_list = ['datalists/woidM36_10k_samples_' strain 'like.txt'];
        filepath = '../../../sworm-model/woidModel/results/paramSampleResults/woids/';
        scaleflag = 'linear';
end
supportRange = [revRate_range; dkdN_range]';

%% Analyse simulations and experiments - or load pre-computed summary statistics

if exist(sumstat_filename,'file')
    load(sumstat_filename)
else
    % analyse simulation data
    [sim_ss_array, sim_file_names, param_return] = f_analyse_sims(sim_file_list,...
        filepath, param_names, num_statistics);
    % Analyse experimental data
    [exp_ss_array, exp_strain_list] = f_analyse_exps({strain},2,num_statistics);
    save(sumstat_filename,'sim_ss_array','sim_file_names','param_return',...
        'exp_ss_array','exp_strain_list','model')
end

%% optimise weightings of summary statistics for model and strain
optimresults_filename = ['optimresults_' model '_' strain '_alpha_' num2str(accept_ratio) '.mat'];
if exist(optimresults_filename,'file')
    load(optimresults_filename)
else
    disp('optimising weights...')
    [weights_optim, min_obj] = f_optim_posterior(exp_ss_array, sim_ss_array,...
        accept_ratio, param_names, param_return, supportRange,scaleflag);
    save(optimresults_filename,'weights_optim','min_obj','strain','model','accept_ratio')
end
%% Obtain distances between each of the experiments and simulations
expsim_dists = f_exp2sim_dist(exp_ss_array, sim_ss_array,weights_optim);

%% Perform parameter inference
[chosen_params, chosen_samples] = f_infer_params(...
    expsim_dists, exp_strain_list,[accept_ratio],{'r_{rev}','dk/d\rho'},param_return,...
    true,supportRange,scaleflag,model);

%% Plot summary statistics of experiments and best samples
exportOptions = struct('Format','eps2','Color','rgb','Width',10,...
    'Resolution',300,'FontMode','fixed','FontSize',10,'LineWidth',1);
for statCtr = 1:2
    sumStatFig = figure;
    for strainCtr = 1:length(exp_strain_list)
        semilogy(exp_ss_array{strainCtr,statCtr+1},'LineWidth',2)
        hold on
    end
    for strainCtr = 1:length(exp_strain_list)
        for ii=1
            semilogy(sim_ss_array{chosen_samples(strainCtr,ii),statCtr+1})
        end
    end
    title(['S_' num2str(statCtr) ', weight ' num2str(weights_optim(statCtr)./sum(weights_optim),2) ],'FontWeight','normal')
    legend([exp_strain_list{1} ' mean'],'best simulations')
    formatAndExportFigure(sumStatFig,['figures/S_' num2str(statCtr) '_' ...
                strain '_alpha_' num2str(accept_ratio) '_' model],exportOptions)
end

%% test coverage
f_test_coverage(chosen_samples,200,ones(size(1./expsim_dists(1,chosen_samples,1))),...
    sim_ss_array,weights_optim,accept_ratio,param_names,param_return,supportRange,true,strain,model)

% %% plot surface of dissimilatirity
% load ../../../sworm-model/woidModel/paramSamples_nSim10000_nParam2.mat
% F = scatteredInterpolant(paramSamples.dkdN,paramSamples.revRateClusterEdge,expsim_dists(:,1));
end