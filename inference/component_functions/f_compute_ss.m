function ss_results = f_compute_ss(in_data, expsim_classifier, ...
    fraction_to_sample, num_statistics)

if nargin<3
    fraction_to_sample = 0.1;
end

% Set up cell array to store function results in
ss_results = cell(1,num_statistics);

ss_results{1} = inf_pcf(in_data, expsim_classifier, fraction_to_sample);

if num_statistics>1
    ss_results{2} = inf_hierarchicalclustering(in_data, expsim_classifier, fraction_to_sample);
    if num_statistics>2
       [ss_results{3}, ss_results{4},ss_results{5}, ss_results{6}] ...
           = inf_positionalmoments(in_data, expsim_classifier, fraction_to_sample); 
        if num_statistics>6
           ss_results{7} = inf_polarorder(in_data, expsim_classifier, fraction_to_sample); 
        end
    end
end
% Extend for as many summary statistics as needed
% e.g. ss_results{n} = ss_function_n(in_data, expsim_classifier)
end