function ss_results = f_compute_ss(in_data, expsim_classifier)
global num_statistics

% Set up cell array to store function results in
ss_results = cell(1,num_statistics-1)

[ss_results{1},ss_results{2},ss_results{3}] = ... 
    inf_clusterdistribution(in_data, expsim_classifier);

ss_results{4} = inf_gr(in_data, expsim_classifier);

% Extend for as many summary statistics as needed
% e.g. ss_results{n} = ss_function_n(in_data, expsim_classifier)
end