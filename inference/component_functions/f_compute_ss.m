function ss_results = f_compute_ss(in_data, expsim_classifier, fraction_to_sample)
global num_statistics

if nargin<3
    fraction_to_sample = 0.1;
end

% Set up cell array to store function results in
ss_results = cell(1,num_statistics-1);

% [ss_results{1},ss_results{2},ss_results{3}] = ... 
%     inf_clusterdistribution(in_data, expsim_classifier);

ss_results{1} = inf_gr(in_data, expsim_classifier, fraction_to_sample);

% Extend for as many summary statistics as needed
% e.g. ss_results{n} = ss_function_n(in_data, expsim_classifier)
end