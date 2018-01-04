function [sim_sstat_array, sim_file_names, param_return] = ...
    f_analyse_sims(sim_file_list, extract_params,num_statistics)

sim_file_names = {};
my_list_file = fopen(sim_file_list);

% Read in the first line of the .txt file
tline = fgetl(my_list_file);
% For every line that still contains text...
while ischar(tline)
    % Store this file location in the initialised array
    sim_file_names{end+1} = tline;
    % Read in the next line
    tline = fgetl(my_list_file);
end
% Close the list file, having read in all of the file names
fclose(my_list_file)

numSims = length(sim_file_names);

% Construct array for storing the outputs of the summary statistics for
% each simulation. Firstly, have to account for the extra column for the
% strain/exp name
sim_sstat_array = cell(numSims, 1+num_statistics);

if iscell(extract_params)
    param_return = zeros(numSims, length(extract_params));
end

filepath = '../../../sworm-model/woidModel/results/paramSampleResults/woidlinos/woidM18paramD2/';
% filepath = '';
% For each simulation file in the list, compute the appropriate summary
% statistics using supplied functions
for simCtr = 1:numSims
    load([filepath sim_file_names{simCtr}]);
    
    % Extract desired parameters if extract_params is a cell array
    if iscell(extract_params)
        for par = 1:length(param_return)
            param_return(simCtr,par) = ...
                eval(strcat('param.', extract_params{par}));
        end
    end
    
    sim_sstat_array{simCtr,1} = sim_file_names{simCtr};
    
    % compute all chosen summary statistics
    fraction_to_sample = min(param.dT*param.saveEvery/3,1); % specifiy fraction of frames to sample
    sstat_results = f_compute_ss(xyarray, 'simulation', fraction_to_sample, num_statistics);
    
    for sstatCtr = 1:length(sstat_results)
        sim_sstat_array{simCtr, sstatCtr+1} = sstat_results{sstatCtr};
    end 
    % display progress
    disp(['analysed ' num2str(100*simCtr/numSims) '% of simulations'])
end
end

