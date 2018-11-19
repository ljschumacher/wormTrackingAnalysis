function [weights_optim, min_obj] = inference_master(model,accept_ratio)
% Inference framework

% issues/to-do:
% - calculation of summary statistics could be sped up by calculating all
% stats within the loop over frames, rather than looping over frames for each stat
addpath('component_functions');

% % set overall parameters
% model = 'log-rods'; % 'rods' or 'worms'
% accept_ratio = 0.005, 0.01, 0.02, or 0.05

switch model
    case 'PRW_4D_wa_r2'
        param_names = {'drdN_rev','dkdN_dwell','dkdN_undwell','f_hapt'};
        num_statistics = 4;
        % load original priors (that have been adjusted to log-scale for parameter 4)
        load(['../../../sworm-model/woidModel/sampling/priors4D_log_M_18_noVolExcl'...
            '_angleNoise_0.05_k_theta_0_slowing_stochastic_bynode_dwell_0.0036_1.1_' ...
            'revdensity_haptotaxis_weighted_additive.mat'],'prior_npr1','supportLimits');
        prior{1} = prior_npr1;
        load(['../../../sworm-model/woidModel/sampling/priors4D_log_M_18_noVolExcl'...
            '_angleNoise_0.0326_k_theta_0_slowing_stochastic_bynode_dwell_0.25_0.45_' ...
            'revdensity_haptotaxis_weighted_additive.mat'],'prior_N2');
        prior{2} = prior_N2;
        % load support Limits from r1 posteriors - the proposal distribution
        load(['inf_results/posteriors_log_PRW_4D_wa_r1_0.1.mat'],'supportLimits')
        % load r1 posteriors - the proposal distribution
        load(['inf_results/posteriors_log_PRW_4D_wa_r1_0.1.mat'],'posterior');
        proposal = posterior;
        sumstat_filename = 'sumstats_PRW_4D_wa_r2.mat';
        sim_file_lists = {'datalists/PRW_4D_wa_r2_27384_samples_npr1like.txt';
            'datalists/PRW_4D_wa_r2_13341_samples_N2like.txt'};
        filepath = {'../../../sworm-model/woidModel/results/woidlinos/paramSamples/PRW_4D_taxis_weighted_additive_r2/npr_1/'; ...
            '../../../sworm-model/woidModel/results/woidlinos/paramSamples/PRW_4D_taxis_weighted_additive_r2/N2/'};
        scaleflag = 'linear';
        supportLimits(1,4)=-4; % fix support limits
    case 'PRW_4D_wa_r1'
        param_names = {'drdN_rev','dkdN_dwell','dkdN_undwell','f_hapt'};
        num_statistics = 4;
        % load old priors
        load(['../../../sworm-model/woidModel/sampling/priors4D_M_18_noVolExcl'...
            '_angleNoise_0.05_k_theta_0_slowing_stochastic_bynode_dwell_0.0036_1.1_' ...
            'revdensity_haptotaxis_weighted_additive_WRONGBANDWIDTH.mat'],'prior_npr1','supportLimits');
        proposal{1} = prior_npr1;
        load(['../../../sworm-model/woidModel/sampling/priors4D_M_18_noVolExcl'...
            '_angleNoise_0.0326_k_theta_0_slowing_stochastic_bynode_dwell_0.25_0.45_' ...
            'revdensity_haptotaxis_weighted_additive_WRONGBANDWIDTH.mat'],'prior_N2');
        proposal{2} = prior_N2;
        % load corrected priors
        load(['../../../sworm-model/woidModel/sampling/priors4D_M_18_noVolExcl'...
            '_angleNoise_0.05_k_theta_0_slowing_stochastic_bynode_dwell_0.0036_1.1_' ...
            'revdensity_haptotaxis_weighted_additive.mat'],'prior_npr1','supportLimits');
        prior{1} = prior_npr1;
        load(['../../../sworm-model/woidModel/sampling/priors4D_M_18_noVolExcl'...
            '_angleNoise_0.0326_k_theta_0_slowing_stochastic_bynode_dwell_0.25_0.45_' ...
            'revdensity_haptotaxis_weighted_additive.mat'],'prior_N2');
        prior{2} = prior_N2;
        
        sumstat_filename = 'sumstats_PRW_4D_wa_r1.mat';
        sim_file_lists = {'datalists/PRW_4D_wa_r1_11214_samples_npr1like.txt';
            'datalists/PRW_4D_wa_r1_1394_samples_N2like.txt'};
        filepath = {'../../../sworm-model/woidModel/results/woidlinos/paramSamples/PRW_4D_taxis_weighted_additive_r1/npr_1/'; ...
            '../../../sworm-model/woidModel/results/woidlinos/paramSamples/PRW_4D_taxis_weighted_additive_r1/N2/'};
        scaleflag = 'linear';
    case 'PRW_2D_pilot'
        param_names = {'drdN_rev','dkdN_dwell'};
        num_statistics = 4;
        sumstat_filename = 'sumstats_PRW_2D_pilot.mat';
        sim_file_lists = {'datalists/PRW_2D_pilot.txt'};
        filepath = {'../../../sworm-model/woidModel/results/woidlinos/floppy/'};
        scaleflag = 'linear';
        supportLimits = [0 0; 1 1];
    case 'worms'
        param_names = {'revRateClusterEdge','dkdN_dwell'};
        num_statistics = 4; % 5th stat, polar order, did not seem to work well
        load('../../../sworm-model/woidModel/paramSamples_wM36_nSim20000_nParam2.mat')
        supportLimits = [revRate_range; dkdN_range]';
        sumstat_filename = ['sumstats_20k_samples_wM36.mat'];
        sim_file_lists = {'datalists/woidM36_20k_samples_npr1like.txt';...
            'datalists/woidM36_20k_samples_N2like.txt'};
        filepath = '../../../sworm-model/woidModel/results/paramSampleResults/woids/';
        scaleflag = 'linear';
end
nParams = length(param_names);
%% Analyse simulations and experiments - or load pre-computed summary statistics

if exist(sumstat_filename,'file')
    load(sumstat_filename)
else
    % analyse simulation data
    [sim_ss_array, sim_file_names, param_return] = f_analyse_sims(sim_file_lists,...
        filepath, param_names, num_statistics);
    save(sumstat_filename,'sim_ss_array','sim_file_names','param_return','model')
    %     % Analyse experimental data
    %     [exp_ss_array, exp_strain_list] = f_analyse_exps({'npr1','N2'},2,num_statistics);
    % load experimental data
    load('sumstats_expmnt.mat')
    save(sumstat_filename,'sim_ss_array','sim_file_names','param_return',...
        'exp_ss_array','exp_strain_list','model')
end
% adjust selected parameters to log-scale
if strcmp(model,'PRW_4D_wa_r2')
    for strainCtr = 1:2
        param_return{strainCtr}(:,4) = log10(param_return{strainCtr}(:,4));
    end
end
%% optimise weightings of summary statistics for model and strain
optimresults_filename = ['optim_results/optimresults_' model '_alpha_' num2str(accept_ratio) '.mat'];
if exist(optimresults_filename,'file')
    load(optimresults_filename)
else
    disp('optimising weights...')
    [weights_optim, min_obj] = f_optim_posterior(exp_ss_array, sim_ss_array,...
        accept_ratio, param_names, param_return, prior,proposal,scaleflag);
    save(optimresults_filename,'weights_optim','min_obj','model','accept_ratio')
end
weights_optim %= [1 1 1 1]
%% Obtain distances between each of the experiments and simulations
expsim_dists = f_exp2sim_dist(exp_ss_array, sim_ss_array,weights_optim);
% check marginals of distances against parameters
for strainCtr = 1:length(sim_ss_array)
    figure
    for paramCtr = 1:nParams
        for statCtr =1:num_statistics
            subplot(num_statistics,nParams,paramCtr + (statCtr - 1)*nParams)
            scatter(param_return{strainCtr}(:,paramCtr),expsim_dists{strainCtr}(:,statCtr+1),'k.')
            hl = refline;
            hl.LineWidth = 2;
            hl.Color = [1 0 0];
            if paramCtr==1, ylabel('distance'), end
            title(['S_' num2str(statCtr)])
            if statCtr==num_statistics, xlabel(param_names{paramCtr}), end
        end
    end
end
%% Perform parameter inference
[chosen_params, chosen_samples] = f_infer_params(expsim_dists,...
    exp_strain_list,[accept_ratio],{'r''','k''_s','k''_f','log_{10}f_t'},...
    param_return,true,supportLimits,scaleflag,model);
% display filenames for best simulations
disp(['best npr-1:' sim_file_names{1}{chosen_samples{1}(1)}])
disp(['best N2:' sim_file_names{2}{chosen_samples{2}(1)}])

%% save posterior distribution
if strcmp(model,'PRW_4D_wa_r1')
    % load original priors (that have been adjusted to log-scale for parameter 4)
    load(['../../../sworm-model/woidModel/sampling/priors4D_log_M_18_noVolExcl'...
        '_angleNoise_0.05_k_theta_0_slowing_stochastic_bynode_dwell_0.0036_1.1_' ...
        'revdensity_haptotaxis_weighted_additive.mat'],'prior_npr1','supportLimits');
    log_prior{1} = prior_npr1;
    load(['../../../sworm-model/woidModel/sampling/priors4D_log_M_18_noVolExcl'...
        '_angleNoise_0.0326_k_theta_0_slowing_stochastic_bynode_dwell_0.25_0.45_' ...
        'revdensity_haptotaxis_weighted_additive.mat'],'prior_N2');
    log_prior{2} = prior_N2;
    chosen_params_log = chosen_params;
    for strainCtr =1:2
        chosen_params_log{strainCtr}(:,4) = log10(chosen_params{strainCtr}(:,4));
        % kernel density estimation
        nSamples = numel(chosen_samples{strainCtr});
        ndims = nParams;
        bw = std(chosen_params_log{strainCtr}).*(4./(ndims + 2)./nSamples).^(1./(ndims + 4)); %Silverman's rule of thumb for the bandwidth
        % impose minimum bandwidth limit - set minimum sigma to be half grid size
        if any(bw==0), error('zero bandwidths'), end
        % % replicate multivariate kernel density estimation using a gaussian mixture model
        % prior = mvksdensity(paramCombis(sampleIndcs,:),paramCombis,'Bandwidth',bw);
        kde_weights = 1./expsim_dists{strainCtr}(chosen_samples{strainCtr},1);
        % adjust weighting for change of priors
        prior_weights = pdf(log_prior{strainCtr},chosen_params_log{strainCtr})...
            ./pdf(proposal{strainCtr},chosen_params{strainCtr});
        posterior{strainCtr} = gmdistribution(chosen_params_log{strainCtr},...
            bw.^2,kde_weights.*prior_weights);
    end
elseif strcmp(model,'PRW_4D_wa_r2')
    for strainCtr =1:2
        % kernel density estimation
        nSamples = numel(chosen_samples{strainCtr});
        ndims = nParams;
        bw = std(chosen_params{strainCtr}).*(4./(ndims + 2)./nSamples).^(1./(ndims + 4)); %Silverman's rule of thumb for the bandwidth
        % impose minimum bandwidth limit - set minimum sigma to be half grid size
        if any(bw==0), error('zero bandwidths'), end
        % % replicate multivariate kernel density estimation using a gaussian mixture model
        % prior = mvksdensity(paramCombis(sampleIndcs,:),paramCombis,'Bandwidth',bw);
        kde_weights = 1./expsim_dists{strainCtr}(chosen_samples{strainCtr},1);
        % adjust weighting for change of priors
        prior_weights = pdf(prior{strainCtr},chosen_params{strainCtr})...
            ./(1e-3 + pdf(proposal{strainCtr},chosen_params{strainCtr})); % use some regularization
        posterior{strainCtr} = gmdistribution(chosen_params{strainCtr},...
            bw.^2,kde_weights.*prior_weights);
    end
end
% save posterior for later sampling
save(['inf_results/posteriors_log_' model '_' num2str(accept_ratio) ...
    '.mat'],'bw','supportLimits','posterior','param_names','num_statistics',...
    'weights_optim','scaleflag')
%% Plot summary statistics of experiments and best samples
exportOptions = struct('Format','eps2','Color','rgb','Width',10,...
    'Resolution',300,'FontMode','fixed','FontSize',10,'LineWidth',1);
plotColors = lines(2);
plotbins = (0.1:0.1:2) - 0.1/2;
for statCtr = 1:2
    sumStatFig = figure;
    for strainCtr = 1:length(sim_ss_array)
        nBins = size(exp_ss_array{strainCtr,statCtr+1},2);
        errorbar(plotbins(1:nBins),mean(exp_ss_array{strainCtr,statCtr+1}),...
            std(exp_ss_array{strainCtr,statCtr+1}),':','LineWidth',2,'Color',plotColors(strainCtr,:))
        hold on
    end
    sumStatFig.Children.YScale = 'log';
    for strainCtr = 1:size(sim_ss_array,1)
        nBins = size(exp_ss_array{strainCtr,statCtr+1},2);
        for ii=1:1
            semilogy(plotbins(1:nBins),sim_ss_array{strainCtr}{chosen_samples{strainCtr}(ii),statCtr+1},'LineWidth',2,'Color',plotColors(strainCtr,:))
        end
    end
    xlabel('r (mm)')
    if statCtr==2
        xlim([0 1.6])
    end
    if statCtr==1
        ylim([0.5 100])
    end
    title(['S_' num2str(statCtr) ])%', weight ' num2str(100*weights_optim(statCtr)./sum(weights_optim),3) '%'])
    legend([exp_strain_list{1} ' mean'],[exp_strain_list{2} ' mean'],[exp_strain_list{1} ' best sim.'],[exp_strain_list{2} ' best sim.'])
    %     legend([exp_strain_list{1} ' mean'],['simulation'])
    formatAndExportFigure(sumStatFig,['figures/S_' num2str(statCtr) ...
        '_alpha_' num2str(accept_ratio) '_' model],exportOptions)
end
for statCtr = 3:4
    sumStatFig = figure;
    for strainCtr = 1:size(sim_ss_array,1)
        subplot(1,2,strainCtr)
        violinplot(cat(1,sim_ss_array{strainCtr}{chosen_samples{strainCtr}(:),statCtr+1}),...
            strainCtr,'ViolinColor',plotColors(strainCtr,:),'ViolinAlpha',1,...
            'BoxColor',plotColors(strainCtr,:),'ShowData',false,'Width',0.1);
        hold on
    end
    for strainCtr = 1:length(exp_strain_list)
        subplot(1,2,strainCtr)
        errorbar(1,mean(exp_ss_array{strainCtr,statCtr+1}),...
            std(exp_ss_array{strainCtr,statCtr+1}),'k+','LineWidth',2)
        set(gca,'YScale','log')
        set(gca,'XTick',[])
        xlabel(exp_strain_list{strainCtr})
    end
    subplot(1,2,1)
    title(['S_' num2str(statCtr) ', weight ' num2str(100*weights_optim(statCtr)./sum(weights_optim),3) '%'])
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
