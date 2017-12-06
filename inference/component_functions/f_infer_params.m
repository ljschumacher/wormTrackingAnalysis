function chosen_params = f_infer_params(expsim_dists, test_params, n_cuts)
global num_statistics
% Set the cutoffs for taking the top n% of simulations e.g to select the
% closest 1% of simulations, use 'n_cuts = [0.01]'. To see the effect that
% using different cutoffs has on the parameter distributions inferred,
% separate values with a comma: 'n_cuts = [0.10,0.05,0.01].

%Create array for storing the parameters of the top n% of simulations
chosen_params = zeros(floor(prod(size(expsim_dists))*max(n_cuts)/num_statistics),...
    length(test_params), length(n_cuts));

% For each of the % cutoffs specified in n_cuts, produce distributions of 
% the parameters
for cutoff = 1:length(n_cuts)
    n = n_cuts(cutoff)
    top_n = floor(prod(size(expsim_dists))*n/num_statistics);
    
    lin = reshape( expsim_dists(:,:,1).' ,1,numel(expsim_dists(:,:,1)));
    A = sort(lin);
    %A = sort(expsim_dists(:,:,1));
    B = (expsim_dists(:,:,1)<=A(top_n)).*expsim_dists(:,:,1);
    
    best_sims = floor(find(B)/length(in_list));
    list_best = {};
    
    for sim = 1:length(best_sims)
        list_best{end+1} = sim_file_names(best_sims(sim)+1);
    end
    
    % Could use the mat file to find the parameters for plotting, 
    % no need to reload each of the simulations. Remains fast without.
    
    for par = 1:length(test_params)
        for i = 1:length(list_best)
            load(list_best{i}{1});
            chosen_params(i,par,cutoff) = eval(strcat('param.', test_params{par}));
        end  
    end    
end

% -------- Producing joint distributions of inferred parameters -------- %
figure;
for cutoff = 1:length(n_cuts)
    to_plot = chosen_params(:,:,cutoff);
    
    % Eliminate redundantrows, where all parameter values are zero
    % Useful when there are multiple, increasingly tight cutoffs
    to_plot = to_plot(any(to_plot~=0,2),:);
    
    subplot(1,length(n_cuts),cutoff)
    [S,AX,BigAx,H,HAx] = plotmatrix(to_plot);
    
    title(['Top ' num2str(n_cuts(cutoff)*100) '% of simulations'])
    
    for i = 1:length(test_params)
        ylabel(AX(i,1),test_params(i))
        xlabel(AX(length(test_params),i),test_params(i))
    end
    
end
end