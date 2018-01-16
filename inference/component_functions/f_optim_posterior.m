function expsim_dists = f_optim_posterior(exp_ss_array, sim_ss_array, exp_strain_list,...
    params, p_cutoff, paramFile)
% Optimise the weights of summary statistics to maximise the Hellinger
% distance between the prior and posterior, and return the the appropriate distances between each of the
% simulations and the experimental references

load(paramFile)
numStrains = length(exp_strain_list);

supportRange = [revRate_range; dkdN_range]';
% construct a kernel density estimate of the prior
[prior, prior_query_points] = ksdensity(paramSample.Variables,...
    'Support',supportRange,'BoundaryCorrection','reflection');
% choose weights

% for given weights, compute the distances
expsim_dists = f_exp2sim_dist(exp_ss_array, sim_ss_array, exp_strain_list, weights);
% for these distances, select fraction of closest parameters
[chosen_params, chosen_samples] = f_infer_params(...
    expsim_dists, exp_strain_list, p_cutoff, paramSamples, false);
for strainCtr = 1:numStrains
    % construct a kernel density estimate of the posterior
    kde_weights = 1./expsim_dists(strainCtr,chosen_samples(strainCtr),1);
    [posti,posti_query_points] = ksdensity(squeeze(chosen_params(strainCtr,:,:)),...
        'Support',supportRange,'BoundaryCorrection','reflection','Weights',kde_weights)
end

end