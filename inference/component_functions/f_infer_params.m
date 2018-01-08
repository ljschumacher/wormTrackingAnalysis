function [chosen_params, chosen_samples] = f_infer_params(expsim_dists,...
    params, p_cutoffs, paramFile, plotResults)
% Set the cutoffs for taking the top p% of simulations e.g to select the
% closest 1% of simulations, use 'p_cutoffs = [0.01]'. To see the effect that
% using different cutoffs has on the parameter distributions inferred,
% separate values with a comma: 'p_cutoffs = [0.04,0.02,0.01].

if nargin<5
    plotResults = false;
end
num_sims = size(expsim_dists,2);
num_strains = size(expsim_dists,1);
load(paramFile)
if plotResults
    load('~/Dropbox/Utilities/colormaps_ascii/increasing_cool/cmap_Blues.txt')
end
%Create array for storing the parameters of the top p% of simulations
chosen_params = zeros(num_strains,floor(num_sims*max(p_cutoffs)),...
    length(params),length(p_cutoffs));
% Create array for returning the original sample indicies of the chosen
% params
chosen_samples = zeros(num_strains,floor(num_sims*max(p_cutoffs)),...
    length(p_cutoffs));

for strainCtr = 1:num_strains
    % For each of the % cutoffs specified in p_cutoffs, produce distributions of
    % the parameters
    disp(['Inferring parameters for strain ' num2str(strainCtr) '/' num2str(num_strains)])
    for cutoffCtr = 1:length(p_cutoffs)
        this_cutoff = p_cutoffs(cutoffCtr);
        num_top_samples = floor(num_sims*this_cutoff);
        
        [sorted_distances, sorted_indeces] = sort(expsim_dists(strainCtr,:,1));
        acceptedSamples_logInd = expsim_dists(strainCtr,:,1)<=sorted_distances(num_top_samples);
        chosen_samples(strainCtr,1:num_top_samples,cutoffCtr) = sorted_indeces(1:num_top_samples);
        
        for paramCtr = 1:length(params)
            thisParamSamples = eval(strcat('paramSamples.', params{paramCtr}));
            chosen_params(strainCtr,1:num_top_samples,paramCtr,cutoffCtr) = ...
                thisParamSamples(acceptedSamples_logInd);
        end
    end
    
    % -------- Producing joint distributions of inferred parameters -------- %
    if plotResults
        figure;
        for cutoffCtr = 1:length(p_cutoffs)
            to_plot = squeeze(chosen_params(strainCtr,:,:,cutoffCtr));
            % Eliminate redundant rows, where all parameter values are zero
            % Occures when there are multiple cutoffs chosen
            to_plot = to_plot(any(to_plot~=0,2),:);  %% unclear if necessary
            
            subplot(1,length(p_cutoffs),cutoffCtr)
            [~,AX,~,~,~] = hplotmatrix(to_plot); %% add KDE, or use hplotmatrix?
            colormap(flipud(cmap_Blues))
            title(['Top ' num2str(p_cutoffs(cutoffCtr)*100) '% of simulations'])
            
            for i = 1:length(params)
                ylabel(AX(i,1),params(i))
                xlabel(AX(length(params),i),params(i))
            end
        end
    end
end
end