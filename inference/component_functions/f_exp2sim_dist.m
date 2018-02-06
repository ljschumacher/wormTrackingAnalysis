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
numSims = size(sim_ss_array,1);
numStrains = size(exp_ss_array,1);
expsim_dists = zeros(numStrains,numSims, 1+num_statistics);

for statCtr = 1:num_statistics
    if isempty(weights)
        % to scale the distances for the scale of the summary
        % statistics, one commonly scales by the standard deviation
        % of each statistic. since we have distributions (ie binned
        % data), we will divide by the deviation for each bin
        normfactor = 1./std(cat(1,sim_ss_array{:,1+statCtr}));
    else
        normfactor = weights(statCtr);%./std(cat(1,sim_ss_array{:,1+statCtr}));
        % for summary statistics that are made up of binned data, we may
        % still want to normalise the scale/range across different bins
        % (so as not to weight bins with higher counts more strongly)
        % or log-transform?
    end
    for strainCtr = 1:numStrains
        exp_data = exp_ss_array{strainCtr,1+statCtr};
        scale_factor = 1;% exp_data;
%         scale_factor(scale_factor==0) = min(scale_factor(scale_factor~=0)); % to normalise by observed summary statistic, take care not to divide by zero
        for simCtr = 1:numSims
            sim_data = sim_ss_array{simCtr,1+statCtr};
            dim_factor = 1./sqrt(length(exp_data)); % correction factor for higher dimensional summary statistics
            % Compute the distance between this simulation and the
            % reference - careful not to take log(0)
            expsim_dists(strainCtr,simCtr,1+statCtr) = norm((log(max(exp_data,eps)) - log(max(sim_data,eps)))./scale_factor...
                .*normfactor).*dim_factor;
            if any(isinf( expsim_dists(strainCtr,simCtr,1+statCtr)))
               1; 
            end
        end
    end
end

% sum the distances from all summary statistics to one combined distance
for strainCtr = 1:numStrains
    for simCtr = 1:numSims
        expsim_dists(strainCtr,simCtr,1) = sum(expsim_dists(strainCtr,simCtr,2:num_statistics+1));
    end
end

end