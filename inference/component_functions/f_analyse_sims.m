function [sim_sstat_array, sim_file_names, param_return] = ...
    f_analyse_sims(sim_file_lists, filepath, param_names,num_statistics,modelclassifier)

if nargin<5 || isempty(modelclassifier)
    modelclassifier = 'simulation';
end

nStrains = length(sim_file_lists);
sim_file_names = cell(nStrains,1);
sim_sstat_array = cell(nStrains,1);
if iscell(param_names)
    param_return = cell(nStrains,1);
end
for strainCtr = 1:nStrains
    this_list_file = fopen(sim_file_lists{strainCtr});
    sim_file_names{strainCtr} = {};
    % Read in the first line of the .txt file
    tline = fgetl(this_list_file);
    % For every line that still contains text...
    while ischar(tline)
        % Store this file location in the initialised array
        sim_file_names{strainCtr}{end+1} = tline;
        % Read in the next line
        tline = fgetl(this_list_file);
    end
    % Close the list file, having read in all of the file names
    fclose(this_list_file)
    
    numSims = length(sim_file_names{strainCtr});
    
    % Construct array for storing the outputs of the summary statistics for
    % each simulation. Firstly, have to account for the extra column for the
    % strain/exp name
    sim_sstat_array{strainCtr} = cell(numSims, 1+num_statistics);
    
    if iscell(param_names)
        param_return{strainCtr} = zeros(numSims, length(param_names));
    end
    
    % For each simulation file in the list, compute the appropriate summary
    % statistics using supplied functions
    for simCtr = 1:numSims
        load([filepath{strainCtr} sim_file_names{strainCtr}{simCtr}]);
        
        % Extract desired parameters if extract_params is a cell array
        if iscell(param_names)
            for parCtr = 1:length(param_names)
                param_return{strainCtr}(simCtr,parCtr) = ...
                    eval(strcat('param.', param_names{parCtr}));
            end
        end
        
        sim_sstat_array{strainCtr}{simCtr,1} = sim_file_names{strainCtr}{simCtr};
        
        % compute all chosen summary statistics
        fraction_to_sample = min(param.dT*param.saveEvery/3,1); % specifiy fraction of frames to sample
        sstat_results = f_compute_ss(xyarray, modelclassifier, fraction_to_sample, num_statistics);
        
        for sstatCtr = 1:length(sstat_results)
            sim_sstat_array{strainCtr}{simCtr, sstatCtr+1} = sstat_results{sstatCtr};
        end
        % display progress
        disp(['analysed ' num2str(100*simCtr/numSims,3) '% of simulations'])
    end
end
end

