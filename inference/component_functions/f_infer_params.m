function [chosen_params, chosen_samples] = f_infer_params(expsim_dists,...
    exp_strain_list, p_cutoffs, param_names, param_values, plotResults, supportLimits,scaleflag,modelstring)
% Set the cutoffs for taking the top p% of simulations e.g to select the
% closest 1% of simulations, use 'p_cutoffs = [0.01]'. To see the effect that
% using different cutoffs has on the parameter distributions inferred,
% separate values with a comma: 'p_cutoffs = [0.04,0.02,0.01].

if nargin<6
    plotResults = false;
end
if nargin<8
    scaleflag = 'linear';
end
if nargin<9
    modelstring = '';
end
num_strains = size(expsim_dists,1);
nParams = length(param_names);
if plotResults
    load('~/Dropbox/Utilities/colormaps_ascii/increasing_cool/cmap_Blues.txt')
    load('~/Dropbox/Utilities/colormaps_ascii/increasing_warm/cmap_Oranges.txt')
    exportOptions = struct('Format','eps2','Color','rgb','Width',14,...
        'Resolution',300,'FontMode','fixed','FontSize',10,'LineWidth',1);
    if strcmp(scaleflag,'log')
        supportLimits = log10(supportLimits);
        for paramCtr=1:nParams
            param_names{paramCtr} = ['log_{10}' param_names{paramCtr}];
        end
    end
end
%Create array for storing the parameters of the top p% of simulations
chosen_params = cell(num_strains,1);
% Create array for returning the original sample indicies of the chosen params
chosen_samples = cell(num_strains,1);

for strainCtr = 1:num_strains
    % For each of the % cutoffs specified in p_cutoffs, produce distributions of
    % the parameters
    num_sims = size(expsim_dists{strainCtr},1);
    chosen_params{strainCtr} = zeros(floor(num_sims*max(p_cutoffs)),nParams,length(p_cutoffs));
    chosen_samples{strainCtr} = zeros(floor(num_sims*max(p_cutoffs)),length(p_cutoffs));
    if ~isempty(exp_strain_list)
        disp(['Inferring parameters for strain ' exp_strain_list{strainCtr}])
    end
    for cutoffCtr = 1:length(p_cutoffs)
        this_cutoff = p_cutoffs(cutoffCtr);
        num_top_samples = floor(num_sims*this_cutoff);
        
        [sorted_distances, sorted_indeces] = sort(expsim_dists{strainCtr}(:,1));
        acceptedSamples_logInd = expsim_dists{strainCtr}(:,1)<=sorted_distances(num_top_samples);
        chosen_samples{strainCtr}(1:num_top_samples,cutoffCtr) = sorted_indeces(1:num_top_samples);
        if plotResults % plot fraction of accepted samples
            distFig = figure;
            H = histogram(sorted_distances);
            hold on
            histogram(sorted_distances(1:num_top_samples),H.BinEdges)
            legend('all samples',[num2str(this_cutoff) ' fraction'])
            xlabel('distance'), ylabel('count'), title(exp_strain_list{strainCtr},'FontWeight','normal')
            formatAndExportFigure(distFig,['figures/diagnostics/distances_' ...
                exp_strain_list{strainCtr} '_alpha_' num2str(this_cutoff) '_' modelstring],...
                exportOptions)
        end
        for paramCtr = 1:nParams
            chosen_params{strainCtr}(1:num_top_samples,paramCtr,cutoffCtr) = ...
                param_values{strainCtr}(acceptedSamples_logInd,paramCtr);
        end
    end
    
    % -------- Producing joint distributions of inferred parameters -------- %
    if plotResults
        for cutoffCtr = 1:length(p_cutoffs)
            postiFig{strainCtr} = figure;
            to_plot = squeeze(chosen_params{strainCtr}(:,:,cutoffCtr));
            % Eliminate redundant rows, where all parameter values are zero
            % Occures when there are multiple cutoffs chosen
            to_plot = to_plot(any(to_plot~=0,2),:);  %% only necessary for multiple cut-offs
            if strcmp(scaleflag,'log')
                to_plot = log10(to_plot);
            end
            subplot(1,length(p_cutoffs),cutoffCtr)
            kde_weights = 1./expsim_dists{strainCtr}(chosen_samples{strainCtr}(:,cutoffCtr),1);
            [h{strainCtr},AX{strainCtr},~,hhist{strainCtr},pax{strainCtr}] ...
                = hplotmatrix(to_plot,[],kde_weights, supportLimits);
            if strainCtr==1
                colormap(flipud(cmap_Blues(1:16:end,:)))
            elseif strainCtr==2
                colormap(flipud(cmap_Oranges(1:16:end,:)))
            end
            title(['Approximate posterior from ' num2str(p_cutoffs(cutoffCtr)*100) '% of simulations'],...
                'FontWeight','normal')
            for paramCtr = 1:nParams
                ylabel(AX{strainCtr}(paramCtr,1),param_names(paramCtr))
                xlabel(AX{strainCtr}(nParams,paramCtr),param_names(paramCtr))
            end
        end
    end
end
if plotResults
    % combine plots from both strains into 1
   plotColors = lines(2);
   for lineCtr = 1:2
       hhist{2}(lineCtr).Color = plotColors(2,:); % make orange
       copyobj(hhist{2}(lineCtr),pax{1}(lineCtr)) % transfer to plot of first strain
   end
   drawnow
   freezeColors(AX{2}(2,1))
   h{1}(2).FaceColor = 'flat'; h{2}(2).FaceColor = 'flat';
   copyobj(h{2}(2),AX{1}(2,1))
   formatAndExportFigure(postiFig{1},['figures/posteriors_' ...
       '_alpha_' num2str(p_cutoffs(cutoffCtr)) '_' modelstring],...
       exportOptions)
   close(postiFig{2})
end
end