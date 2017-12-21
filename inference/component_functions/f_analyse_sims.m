function [sim_sstat_array, sim_file_names, param_return] = ...
    analyse_sims(sim_file_list, extract_params)
global num_statistics

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

% Construct array for storing the outputs of the summary statistics for
% each simulation. Firstly, have to account for the extra column for the
% strain/exp name
sim_sstat_array = cell(length(sim_file_names), num_statistics);

if iscell(extract_params)
    param_return = zeros(length(sim_file_names), length(extract_params));
end

% For each simulation file in the list, compute the appropriate summary
% statistics using supplied functions
for sim = 1:length(sim_file_names)
    progress = sim/length(sim_file_names)
    load(sim_file_names{sim});
    
    % Extract desired parameters if extract_params is a cell array
    if iscell(extract_params)
        for par = 1:length(param_return)
            param_return(sim,par) = ...
                eval(strcat('param.', extract_params{par}));
        end
    end
    
    sim_sstat_array{sim,1} = sim_file_names{sim};
    % Read in the appropriate data
    data = xyarray;
    
    % Then compute all chosen summary statistics
    sstat_results = f_compute_ss(in_data, 'simulation')
    
    for sstatCtr = 1:length(sstat_results)
        sim_sstat_array{sim, sstatCtr+1} = sstat_results{sstatCtr}
    end 
    
end
end

