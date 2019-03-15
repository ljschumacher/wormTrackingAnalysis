function [weights_optim,min_obj] = f_optim_posterior(exp_ss_array, sim_ss_array,...
    p_cutoff, param_names, param_values, prior,proposal,scaleflag,dist_method)
% Optimise the weights of summary statistics to maximise the Hellinger
% distance between the prior and posterior, and return the the appropriate distances
% between each of the simulations and the experimental references
% issues / to-do:
% - we may be able to improve performance by not evaluating the prior again
% for each strain, if we can rely in ksdensity to always query the same points
if nargin<7
    scaleflag = 'linear';
end
if strcmp(scaleflag,'log')
    for strainCtr = 1:length(param_values)
        param_values{strainCtr} = log10(param_values{strainCtr});
    end
end
numStats = size(exp_ss_array,2) - 1;

% set up objecive function distance
lambda = 1e-4; % regularization parameter
% define the objective function for joint optimisation over both strains
L = @(x) hellinger(x,exp_ss_array,sim_ss_array,p_cutoff,param_names,...
    param_values,prior,proposal,dist_method) ...
    + lambda*norm(x,1); % include regularisation
nIter = max(10*numStats,50);
initial_weights = lhsdesign(nIter,numStats);% just specify the extreme values in here, the rest should be filled in by ga randomly
initial_weights(initial_weights==0) = 1e-14; % to satisfy the optimisation bounds

%% use the ga global optimization toolbox solver
options = optimoptions('ga','PopulationSize',nIter,'InitialPopulationMatrix',initial_weights,'UseParallel',true,'Display','iter',...
    'HybridFcn',@fmincon); % 'UseParallel', true
[weights_optim,min_obj] = ga(L,numStats,[],[],[],[],eps*ones(numStats,1),100*ones(numStats,1),[],options);

% normalise weights (for convenience only)
weights_optim = weights_optim./sum(weights_optim);

% %% use the particle swarm optimizer
% options = optimoptions('particleswarm','InitialSwarmMatrix',initial_weights,'Display','iter',...
%     'HybridFcn',@fmincon);
% [wsw,Lsw] = particleswarm(L,numStats,1e-2*ones(numStats,1),10*ones(numStats,1),options);

% %% use simulated annealing
% options = optimoptions('simulannealbnd','Display','iter');
% [wsa,Lsa] = simulannealbnd(L,mean(initial_weights),zeros(numStats,1),[],options);
end

function L = hellinger(weights,exp_ss_array,sim_ss_array,p_cutoff,...
    param_names,param_values,prior,proposal,dist_method)
% for given weights, compute the distances
if strcmp(dist_method,'log')
    expsim_dists = f_exp2sim_dist(exp_ss_array, sim_ss_array, weights);
elseif strcmp(dist_method,'log_0inStd')
    expsim_dists = f_exp2sim_dist_0inStd(exp_ss_array, sim_ss_array, weights);
end
% for these distances, select fraction of closest parameters
[chosen_params, chosen_samples] = f_infer_params(...
    expsim_dists, [], p_cutoff, param_names, param_values, false, []);
% construct a kernel density estimate of the posterior
nStrains = size(sim_ss_array,1);
nParams = length(param_names);
H = 0;
for strainCtr = 1:nStrains
    kde_weights = 1./expsim_dists{strainCtr}(chosen_samples{strainCtr},1);
    % construct a gmmodel (because ksdensity doesn't work with more
    % than two parameter dimensions)
    bandWidth = std(chosen_params{strainCtr}).*(4./(nParams + 2)...
        ./size(chosen_params{strainCtr},1)).^(1./(nParams + 4)); %Silverman's rule of thumb for the bandwidth
    % adjust weighting for change of sampling distribution, if applicable
    prior_weights = pdf(prior{strainCtr},chosen_params{strainCtr})...
        ./pdf(proposal{strainCtr},chosen_params{strainCtr});
    posti = gmdistribution(chosen_params{strainCtr},bandWidth.^2,kde_weights.*prior_weights);
    % calculate the hellinger distance between the prior and posterior distributions
    H = H + 1./sqrt(2).*norm(sqrt(pdf(posti,param_values{strainCtr}))...
        - sqrt(pdf(prior{strainCtr},param_values{strainCtr})));
end
% give back the negative of the distance, as this is the objective to
% minimize
L = -H;
end

function [c, ceq] = weights_constraint(weights) % may not need this
c = [];
ceq = sum(weights,2) -1; % vectorized
end
