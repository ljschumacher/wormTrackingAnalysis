function [wga,Lga,wsw,Lsw] = f_optim_posterior(exp_ss_array, sim_ss_array,...
    p_cutoff, paramFile)
% Optimise the weights of summary statistics to maximise the Hellinger
% distance between the prior and posterior, and return the the appropriate distances between each of the
% simulations and the experimental references

load(paramFile)
numStats = size(sim_ss_array,2) - 1;

supportRange = [revRate_range; dkdN_range]';
% construct a kernel density estimate of the prior
[prior, prior_query_points] = ksdensity(paramSamples.Variables,...
    'Support',supportRange,'BoundaryCorrection','reflection');
% set up objecive function distance
lambda = 1e-4; % regularization parameter
L = @(x) hellinger(x,prior,prior_query_points,...
    exp_ss_array,sim_ss_array,p_cutoff,paramSamples,supportRange) ...
    + lambda*norm(x,1); % include regularisation
nIter = 10*numStats;
initial_weights = rand(nIter,numStats);

%% use the ga global optimization toolbox solver
options = optimoptions('ga','InitialPopulationMatrix',initial_weights,'Display','iter',...
    'HybridFcn',@fmincon);
[wga,Lga] = ga(L,numStats,[],[],[],[],zeros(numStats,1),[],[],options);

%% use the particle swarm optimizer
options = optimoptions('particleswarm','InitialSwarmMatrix',initial_weights,'Display','iter',...
    'HybridFcn',@fmincon);
[wsw,Lsw] = particleswarm(L,numStats,1e-2*ones(numStats,1),10*ones(numStats,1),options);

% %% use simulated annealing
% options = optimoptions('simulannealbnd','Display','iter');
% [wsa,Lsa] = simulannealbnd(L,mean(initial_weights),zeros(numStats,1),[],options);
end

function L = hellinger(weights,prior,prior_query_points,...
    exp_ss_array,sim_ss_array,p_cutoff,paramSamples,supportRange)
% for given weights, compute the distances
expsim_dists = f_exp2sim_dist(exp_ss_array, sim_ss_array, weights);
% for these distances, select fraction of closest parameters
[chosen_params, chosen_samples] = f_infer_params(...
    expsim_dists, [], p_cutoff, paramSamples, false);
% construct a kernel density estimate of the posterior
kde_weights = 1./expsim_dists(1,chosen_samples,1);
[posti,posti_query_points] = ksdensity(squeeze(chosen_params(1,:,:)),...
    'Support',supportRange,'BoundaryCorrection','reflection','Weights',kde_weights);
% check that density estimates are evaluated at the same points
assert(all(all(prior_query_points==posti_query_points)))
% calculate the hellinger distance between the prior and posterior distributions
H = 1./sqrt(2).*norm(sqrt(posti) - sqrt(prior));
% give back the negative of the distance, as this is the objective to
% minimize
L = -H;
end

function [c, ceq] = weights_constraint(weights) % may not need this
c = [];
ceq = sum(weights) -1;
end