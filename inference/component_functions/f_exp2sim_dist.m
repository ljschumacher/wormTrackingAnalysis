function expsim_dists = f_exp2sim_dist(exp_ss_array, sim_ss_array, in_list)
global num_statistics

% Then, compute the appropriate distances between each of the
% simulations and the experimental references

% Constuct a n-by-m-by-l+1 matrix for containing the distances between each 
% experiment (n) and each simulation (m) for each of the summary statistics
% computed (l)

expsim_dists = zeros(length(in_list),size(sim_ss_array,1), num_statistics);

for i = 1:length(in_list)
    for j = 1:size(sim_ss_array,1)
        
        if i ==1 & j ==1
            running_dists = zeros(1,size(exp_ss_array,2)-1,length(in_list));
        end
        
        exp_data = exp_ss_array(i,2:end);
        sim_data = sim_ss_array(j,2:end);
        
        for ss = 1:length(sim_data)
            if length(exp_data{ss})>1
                exp_data{ss} = norm(exp_data{ss}-sim_data{ss});
                sim_data{ss} = 0;
            end
        end
        
        % Compute the distance between this simulation and the reference
        dists = cell2mat(exp_data)-cell2mat(sim_data);
        running_dists(:,:,i) = running_dists(:,:,i)+dists;
        expsim_dists(i,j,:) = horzcat(norm(dists),dists);
    end
end

% Normalization: divide the distances from each ss by the mean for that
% strain. This ensures that each summary statistic contributes is weighted
% similarly in the distance function. Else, summary statistics at a larger
% scale i.e. % of worms in cluster, would dwarf others when the mean
% euclidean distance is calculated.

for strain = 1:length(in_list)
    for ss = 2:num_statistics
        strain_dist_means(strain, ss-1) = mean(expsim_dists(strain,:,ss));
        expsim_dists(strain,:,ss) = expsim_dists(strain,:,ss)...
            ./mean(expsim_dists(strain,:,ss));  
    end
    
    for sim = 1:size(sim_ss_array,1)
        expsim_dists(strain,sim,1) = sum(expsim_dists(strain,sim,2:num_statistics));
    end
end


% Consider the distance composition to check whether this normalization of
% the summary statistic weights has been successful.
sim_strain_vars = zeros(length(in_list),num_statistics-1);
strain_dist_means = zeros(length(in_list),num_statistics-1);

for strain = 1:length(in_list)
    for ss = 2:num_statistics
        sim_strain_vars(strain, ss-1) = var(expsim_dists(strain,:,ss));
    end
end

% To ensure pie is full
if sum(sim_strain_vars)<=1
    sim_strain_vars(:) = sim_strain_vars(:).*(1/sum(sim_strain_vars));
end

pie_data = zeros(strain,num_statistics-1);
figure;
for strain = 1:length(in_list)
    for ss = 1:num_statistics-1
        pie_data(strain,ss) = sum(expsim_dists(strain,:,ss+1));
    end
        
    subplot(2,length(in_list)+1,strain)
    pie(pie_data(strain,:))
    title(in_list(strain), 'interpreter','none')
    
    subplot(2,length(in_list)+1,strain+length(in_list)+1)
    pie(sim_strain_vars(strain,:))
    title('var in distances')
end
subplot(2,length(in_list)+1,strain+1)
pie(sum(pie_data))
title('Average across strains')
end