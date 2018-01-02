function chosen_params = f_infer_params(expsim_dists, params, p_cutoffs)
global num_statistics
% Set the cutoffs for taking the top p% of simulations e.g to select the
% closest 1% of simulations, use 'p_cutoffs = [0.01]'. To see the effect that
% using different cutoffs has on the parameter distributions inferred,
% separate values with a comma: 'p_cutoffs = [0.04,0.02,0.01].

%Create array for storing the parameters of the top p% of simulations
chosen_params = zeros(floor(prod(size(expsim_dists))*max(p_cutoffs)/num_statistics),...
    length(params), length(p_cutoffs));

% For each of the % cutoffs specified in p_cutoffs, produce distributions of 
% the parameters
for cutoffCtr = 1:length(p_cutoffs)
    this_cutoff = p_cutoffs(cutoffCtr)
    top_samples = floor(prod(size(expsim_dists))*this_cutoff/num_statistics);
            %%% THIS NEEDS TO CHECKED IN DEBUG MODE

    lin = reshape( expsim_dists(:,:,1).' ,1,numel(expsim_dists(:,:,1)));
    A = sort(lin);
    %A = sort(expsim_dists(:,:,1));
    B = (expsim_dists(:,:,1)<=A(top_samples)).*expsim_dists(:,:,1);
    
    best_sims = floor(find(B)/length(in_list));
    list_best = {};
    
    for sim = 1:length(best_sims)
        list_best{end+1} = sim_file_names(best_sims(sim)+1);
    end
    
    % Could use the mat file to find the parameters for plotting, 
    % no need to reload each of the simulations. Remains fast without.
    
    for paramCtr = 1:length(params)
        for i = 1:length(list_best)
            load(list_best{i}{1});
            chosen_params(i,paramCtr,cutoffCtr) = eval(strcat('param.', params{paramCtr}));
        end  
    end    
    %%%
end

% -------- Producing joint distributions of inferred parameters -------- %
figure;
for cutoffCtr = 1:length(p_cutoffs)
    to_plot = chosen_params(:,:,cutoffCtr);
                %%% THIS NEEDS TO CHECKED IN DEBUG MODE

    % Eliminate redundantrows, where all parameter values are zero
    % Useful when there are multiple, increasingly tight cutoffs
    to_plot = to_plot(any(to_plot~=0,2),:);
    
    subplot(1,length(p_cutoffs),cutoffCtr)
    [S,AX,BigAx,H,HAx] = plotmatrix(to_plot); %% add KDE?
    
    title(['Top ' num2str(p_cutoffs(cutoffCtr)*100) '% of simulations'])
    
    for i = 1:length(params)
        ylabel(AX(i,1),params(i))
        xlabel(AX(length(params),i),params(i))
    end
    %%%
end
end