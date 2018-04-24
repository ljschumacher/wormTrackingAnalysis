function  projected_complex_dists = apply_dist_map( ...
    expfullsimple_dists, full_sim_params, smooth_map, smooth_p1, smooth_p2)

projected_complex_dists = expfullsimple_dists(:,:,1)

for i = 1:length(projected_complex_dists(:,:,1))
    
    %Recover the value of revRate and vs
    revRate_value = full_sim_params(i,1);
    vs_value = full_sim_params(i,2);
    
    [M,revRate_map_index] = min(abs(smooth_p1(:)-revRate_value));
    [M,vs_map_index] = min(abs(smooth_p2(:)-vs_value));
    
    dist_scaler = smooth_map(vs_map_index, revRate_map_index);
    
    projected_complex_dists(1,i,1) = projected_complex_dists(1,i,1) * dist_scaler;
    
end
end