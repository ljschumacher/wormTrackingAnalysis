function [smooth_map,smooth_p1_used,smooth_p2_used] = f_create_map(...
    simple_param_return, complex_param_return, ... 
    expsimple_dists, expcomplex_dists, mapping_parameters)
 %% THIS NEEDS EDITING
list_parameter_values_used = zeros(size(simple_param_return,1),length(mapping_parameters),2);

for mapping_point = 1:size(complex_param_return,1)
    
    for parCtr = 1:length(mapping_parameters)
        list_parameter_values_used(mapping_point,parCtr,1) = ...
            simple_param_return(mapping_point,parCtr);
    end
    
    for parCtr = 1:length(mapping_parameters)
        list_parameter_values_used(mapping_point,parCtr,2) = ...
            complex_param_return(mapping_point,parCtr);
    end
end

% Check to see if the parameters for the complex and simple simulations
% are identical, from the order inputted when files were read
if unique(list_parameter_values_used(:,:,1) == list_parameter_values_used(:,:,2)) ~= 1
    report_error = 'Simple and complex simulations have not been inputted in the same order'
else
    
    % Obtain the ratios between the simple and experimental distances
    rough_scales = expcomplex_dists(:,:,1) ./ expsimple_dists(:,:,1);
    
    % Sort the parameter pairs first by vs, then by revRate
    parameter_pairs_raw = list_parameter_values_used(:,:,1);
    [parameter_pairs_sorted, sorting_index] = sortrows(parameter_pairs_raw, [1,2]);
    
    % Order to scale factors in line with the parameters
    rough_scales_ordered = rough_scales(sorting_index);
    
    % Convert linear, ordered list to matrix
    rough_map = flipud(reshape(rough_scales_ordered, ...
        [length(unique(parameter_pairs_sorted(:,2))) length(unique(parameter_pairs_sorted(:,1)))]));
    
    % Check that the sorting and reshaping corretly orders the parameter
    % space
    parameter_pairs_revRate = parameter_pairs_raw(:,1);
    param_map_revRate = flipud(reshape(parameter_pairs_revRate(sorting_index), ...
        [length(unique(parameter_pairs_sorted(:,2))) length(unique(parameter_pairs_sorted(:,1)))]));
    
    parameter_pairs_vs = parameter_pairs_raw(:,2);
    param_map_vs = flipud(reshape(parameter_pairs_vs(sorting_index), ...
        [length(unique(parameter_pairs_sorted(:,2))) length(unique(parameter_pairs_sorted(:,1)))]));
    
    % How many intervals are wanted between cells in the rough map
    num_smoothing_intervals = 2;
    
    %Smooth the rough map with interpolation
    smooth_map = interp2(rough_map, num_smoothing_intervals);
    
    % Also smooth the corresponding parameter maps
    smooth_param_map_revRate = interp2(param_map_revRate, num_smoothing_intervals);
    smooth_param_map_vs = interp2(param_map_vs, num_smoothing_intervals);
    
    smooth_p1_used = unique(smooth_param_map_revRate);
    smooth_p2_used = unique(smooth_param_map_vs);
end
end
