function [weights_optim, min_obj] = inference_master(model,accept_ratio)
% Inference framework

% issues/to-do:
% - calculation of summary statistics could be sped up by calculating all
% stats within the loop over frames, rather than looping over frames for each stat
addpath('component_functions');

% % set overall parameters
% model = 'log-rods'; % 'rods' or 'worms'
% accept_ratio = 0.005, 0.01, 0.02, or 0.05
param_names = {'revRateClusterEdge','dkdN_dwell'};

switch model
    case 'rods'
        num_statistics = 4;
        load('../../sworm-model/woidModel/paramSamples_nSim20000_nParam2.mat')
        sumstat_filename = ['sumstats_20ksamples_wlM18.mat'];
        sim_file_lists = {'datalists/woidM18_20k_samples_npr1like.txt';...
                          'datalists/woidM18_20k_samples_N2like.txt'};
        filepath = '../../../sworm-model/woidModel/results/paramSampleResults/woidlinos/woidM18paramD2/';
        scaleflag = 'linear';
    case 'log-rods'
        num_statistics = 4;
        load('../../sworm-model/woidModel/paramSamples_log_nSim50000_nParam2.mat')
        sumstat_filename = ['sumstats_50klogsamples_wlM18.mat'];
        sim_file_lists = {'datalists/woidM18_50k_logsamples_npr1like.txt';...
                          'datalists/woidM18_50k_logsamples_N2like.txt'};
        filepath = '../../../sworm-model/woidModel/results/paramSampleResults/paramSamplesLog/woidlinos/';
        scaleflag = 'log';
    case 'worms'
        num_statistics = 4; % 5th stat, polar order, did not seem to work well
        load('../../sworm-model/woidModel/paramSamples_wM36_nSim20000_nParam2.mat')
        sumstat_filename = ['sumstats_20k_samples_wM36.mat'];
        sim_file_lists = {'datalists/woidM36_20k_samples_npr1like.txt';...
                         'datalists/woidM36_20k_samples_N2like.txt'};
        filepath = '../../../sworm-model/woidModel/results/paramSampleResults/woids/';
        scaleflag = 'linear';
end
supportRange = [revRate_range; dkdN_range]';

%% Analyse simulations and experiments - or load pre-computed summary statistics

if exist(sumstat_filename,'file')
    load(sumstat_filename)
else
    % analyse simulation data
    [sim_ss_array, sim_file_names, param_return] = f_analyse_sims(sim_file_lists,...
        filepath, param_names, num_statistics);
    % Analyse experimental data
    [exp_ss_array, exp_strain_list] = f_analyse_exps({'npr1','N2'},2,num_statistics);
    save(sumstat_filename,'sim_ss_array','sim_file_names','param_return',...
        'exp_ss_array','exp_strain_list','model')
end

%% optimise weightings of summary statistics for model and strain
optimresults_filename = ['optim_results/optimresults_' model '_alpha_' num2str(accept_ratio) '.mat'];
if exist(optimresults_filename,'file')
    load(optimresults_filename)
else
    disp('optimising weights...')
    [weights_optim, min_obj] = f_optim_posterior(exp_ss_array, sim_ss_array,...
        accept_ratio, param_names, param_return, supportRange,scaleflag);
    save(optimresults_filename,'weights_optim','min_obj','model','accept_ratio')
end
weights_optim
%% Obtain distances between each of the experiments and simulations
expsim_dists = f_exp2sim_dist(exp_ss_array, sim_ss_array,weights_optim);

%% Perform parameter inference
[chosen_params, chosen_samples] = f_infer_params(...
    expsim_dists, exp_strain_list,[accept_ratio],{'r_{rev}','dk/d\rho'},param_return,...
    true,supportRange,scaleflag,model);

%% Plot summary statistics of experiments and best samples
exportOptions = struct('Format','eps2','Color','rgb','Width',10,...
    'Resolution',300,'FontMode','fixed','FontSize',10,'LineWidth',1);
plotColors = lines(2);
plotbins = (0.1:0.1:2) - 0.1/2;
for statCtr = 1:2
    sumStatFig = figure;
    for strainCtr = 1:length(exp_strain_list)
        errorbar(plotbins,mean(exp_ss_array{strainCtr,statCtr+1}),...
            std(exp_ss_array{strainCtr,statCtr+1}),':','LineWidth',2,'Color',plotColors(strainCtr,:))
        hold on
    end
    sumStatFig.Children.YScale = 'log';
    for strainCtr = 1:length(exp_strain_list)
        for ii=1
            semilogy(plotbins,sim_ss_array{strainCtr}{chosen_samples{strainCtr}(ii),statCtr+1},'LineWidth',2,'Color',plotColors(strainCtr,:))
        end
    end
    xlabel('r (mm)')
    title(['S_' num2str(statCtr) ', weight ' num2str(weights_optim(statCtr)./sum(weights_optim),2) ],'FontWeight','normal')
    legend([exp_strain_list{1} ' mean'],[exp_strain_list{2} ' mean'],[exp_strain_list{1} ' best simulation'],[exp_strain_list{2} 'best simulation'])
    formatAndExportFigure(sumStatFig,['figures/S_' num2str(statCtr) ...
        '_alpha_' num2str(accept_ratio) '_' model],exportOptions)
end

% make table or so of summary stat weightings?
%% test coverage
% f_test_coverage(chosen_samples,200,ones(size(1./expsim_dists(1,chosen_samples,1))),...
%     sim_ss_array,weights_optim,accept_ratio,param_names,param_return,supportRange,true,strain,model)

% %% plot surface of dissimilatirity
% figure
% xq = logspace(-1,1,100);
% yq = logspace(-3,0,100);
% [XQ, YQ] = meshgrid(xq,yq);
% for distCtr = 1:(num_statistics+1)
%     subplot(1,num_statistics+1,distCtr)
%     F = RegularizeData3D(paramSamples.revRateClusterEdge,paramSamples.dkdN,squeeze(expsim_dists(1,:,distCtr))',...
%     xq,yq,'smoothness',2e-2,'interp','bicubic','overlap',0.2);
%     contourf(xq,yq,F,100,'EdgeColor','none')
%     set(gca,'XScale','log','YScale','log')
%     if distCtr>1
%         title(['S_' num2str(distCtr-1) ', w=' num2str(weights_optim(distCtr-1),2)], 'FontWeight', 'normal')
%     end
% end
end
