function [] = f_test_coverage(chosen_samples,num_test_samples,...
    sim_ss_array,sum_stat_weights,accept_ratio,paramSamples,plot_flag)
% test coverage to assess quality of posterios
% (see Prangle 2013, van der Vaart 2017)
numParams = size(paramSamples,2);
test_samples = randsample(chosen_samples,num_test_samples,'false');
proportion_under = NaN(num_test_samples,numParams);
for sampleCtr = 1:num_test_samples
    % calculate distances to remaining samples
    this_ss_array = sim_ss_array(test_samples(sampleCtr),:);
    this_dists = f_exp2sim_dist(this_ss_array,sim_ss_array,sum_stat_weights);
    % exclude self from close samples
    this_dists(1,test_samples(sampleCtr),1) = Inf;
    % accept closest simulations
    [this_chosen_params, ~] = f_infer_params(this_dists, [], accept_ratio,...
        paramSamples, false);
    % calculate proportions
    this_true_params = paramSamples{test_samples(sampleCtr),:};
    for paramCtr = 1:numParams
        proportion_under(sampleCtr,paramCtr) = ...
            mean(this_chosen_params(1,:,paramCtr)<=this_true_params(paramCtr));
    end
end

% plot results
if plot_flag
    figure
    for paramCtr = 1:numParams
        subplot(1,numParams,paramCtr)
        histogram(proportion_under(:,paramCtr),10,...
            'Normalization','Probability','DisplayStyle','stairs')
        % test difference from uniform using ks-test
        [~, p] = kstest(proportion_under(:,paramCtr),'CDF',makedist('Uniform'));
        title(['p_{ks} = ' num2str(p,2)])
    end
end

end

