function [] = f_test_coverage(chosen_samples,num_test_samples,sample_weights,...
    sim_ss_array,sum_stat_weights,accept_ratio,param_names,param_values,supportLimits,...
    plot_flag,strain,modelstring)
% test coverage to assess quality of posteriors
% (see Prangle 2013, van der Vaart 2017)
numParams = length(param_names);
test_samples = randsample(chosen_samples,num_test_samples,'true',sample_weights);
proportion_under = NaN(num_test_samples,numParams);
if plot_flag
    exportOptions = struct('Format','eps2','Color','rgb','Width',10,...
    'Resolution',300,'FontMode','fixed','FontSize',10,'LineWidth',1);
end
for sampleCtr = 1:num_test_samples
    % calculate distances to remaining samples
    this_ss_array = sim_ss_array(test_samples(sampleCtr),:);
    this_dists = f_exp2sim_dist(this_ss_array,sim_ss_array,sum_stat_weights);
    % exclude self from close samples
    this_dists(1,test_samples(sampleCtr),1) = Inf;
    % accept closest simulations
    [this_chosen_params, ~] = f_infer_params(this_dists, [], accept_ratio,...
         param_names, param_values, false, supportLimits);
    % calculate proportions
    this_true_params = param_values(test_samples(sampleCtr),:);
    for paramCtr = 1:numParams
        proportion_under(sampleCtr,paramCtr) = ...
            mean(this_chosen_params(1,:,paramCtr)<=this_true_params(paramCtr));
    end
end

% plot results
if plot_flag
    covFig = figure;
    for paramCtr = 1:numParams
        subplot(1,numParams,paramCtr)
        histogram(proportion_under(:,paramCtr),20,...
            'Normalization','Probability','DisplayStyle','stairs')
        % test difference from uniform using ks-test
        [~, p] = kstest(proportion_under(:,paramCtr),'CDF',makedist('Uniform'));
        title([param_names{paramCtr} ', p_{ks} = ' num2str(p,2)],'FontWeight','normal')
    end
    formatAndExportFigure(covFig,['figures/diagnostics/coverage_', ...
                strain '_alpha_' num2str(accept_ratio) '_' modelstring],exportOptions)
end

end

