function expsim_dists = f_exp2sim_dist(exp_ss_array, sim_ss_array, weights)
% Compute the appropriate distances between each of the
% simulations and the experimental references

% Constuct a n-by-m-by-l+1 matrix for containing the distances between each
% experiment (n) and each simulation (m) for each of the summary statistics
% computed (l)

num_statistics = size(exp_ss_array,2)-1;

if nargin<3
    weights=[];
end
numStrains = size(sim_ss_array,1);
expsim_dists = cell(numStrains,1);

for strainCtr = 1:numStrains
    numSims = size(sim_ss_array{strainCtr},1);
    expsim_dists{strainCtr} = zeros(numSims, 1+num_statistics);
    for statCtr = 1:num_statistics
        exp_data = exp_ss_array{strainCtr,1+statCtr};
        assert(~any(exp_data(:)==0),'zero-experimental data, cannot compute log-distances')
        exp_err_mean = std(exp_data)./(sqrt(size(exp_data,1)));
        for simCtr = 1:numSims
            sim_data = sim_ss_array{strainCtr}{simCtr,1+statCtr};
            dim_factor = 1./sqrt(size(exp_data,2)); % correction factor for higher dimensional summary statistics
            % Compute the distance between this simulation and the
            % reference - careful not to take log(0)
%             sim_data(sim_data==0) = exp_err_mean(sim_data==0); % to prevent log(0)
            expsim_dists{strainCtr}(simCtr,1+statCtr) = sum(vecnorm(...
                (log(exp_data) - log(sim_data))... % take scaled difference of all observed values of this summary stat and this simulated one
                ,2,2)... % take norm for each expmntl sample
                .*dim_factor... % correct for dim of summary stat
                );% sum this distance over expmntl samples
            % add the distance to the total from all summary statistics
            expsim_dists{strainCtr}(simCtr,1) = expsim_dists{strainCtr}(simCtr,1)...
                + weights(statCtr).*expsim_dists{strainCtr}(simCtr,1+statCtr); % weight summary statistic
            if statCtr==num_statistics&&expsim_dists{strainCtr}(simCtr,1)==0
                warning('zero distance btw expmnt and simulation');
            end
        end
    end
end

end