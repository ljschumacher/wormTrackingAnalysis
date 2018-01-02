function chosen_params = f_infer_params(expsim_dists, params, p_cutoffs, paramFile)
% Set the cutoffs for taking the top p% of simulations e.g to select the
% closest 1% of simulations, use 'p_cutoffs = [0.01]'. To see the effect that
% using different cutoffs has on the parameter distributions inferred,
% separate values with a comma: 'p_cutoffs = [0.04,0.02,0.01].

num_statistics = size(expsim_dists,3);
load(paramFile)

%Create array for storing the parameters of the top p% of simulations
chosen_params = zeros(floor(prod(size(expsim_dists))*max(p_cutoffs)/num_statistics),...
    length(params), length(p_cutoffs));

% For each of the % cutoffs specified in p_cutoffs, produce distributions of 
% the parameters
for cutoffCtr = 1:length(p_cutoffs)
    this_cutoff = p_cutoffs(cutoffCtr)
    num_top_samples = floor(prod(size(expsim_dists))*this_cutoff/num_statistics);

    sorted_distances = sort(expsim_dists(:,:,1)); % this syntax will sort distances across strains, which may not be what we want if we call the function with multiple strains
    acceptedSamples_logInd = expsim_dists(:,:,1)<=sorted_distances(num_top_samples);
    
    for paramCtr = 1:length(params)
        thisParamSamples = eval(strcat('paramSamples.', params{paramCtr}));
        chosen_params(:,paramCtr,cutoffCtr) = thisParamSamples(acceptedSamples_logInd);
    end    
end

% -------- Producing joint distributions of inferred parameters -------- %
figure;
for cutoffCtr = 1:length(p_cutoffs)
    to_plot = chosen_params(:,:,cutoffCtr);

    % Eliminate redundantrows, where all parameter values are zero
    % Useful when there are multiple, increasingly tight cutoffs
    to_plot = to_plot(any(to_plot~=0,2),:);  %% unclear if necessary
    
    subplot(1,length(p_cutoffs),cutoffCtr)
    [~,AX,~,~,~] = plotmatrix(to_plot); %% add KDE, or use hplotmatrix?
    
    title(['Top ' num2str(p_cutoffs(cutoffCtr)*100) '% of simulations'])
    
    for i = 1:length(params)
        ylabel(AX(i,1),params(i))
        xlabel(AX(length(params),i),params(i))
    end
end
end