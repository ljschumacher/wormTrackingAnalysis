%% Function to output distribution of cluster size from many frames
function [wpc_hist, branch_hist, perc_worms_in_clust] ...
    = inf_clusterdistribution(data, format_in)

% Specify the number of frames to sample
num_samples = 500;
pop_thresh = 2;
cluster_width_cutoffs = [2];

wpc_counts = zeros(1, num_samples);
per_worms_in_clust = zeros(1, num_samples);

if format_in == 'simulation'
    
    % Get the dimensions of the dataframe
    dims = size(data);
    data = data(:,1,:,:);
    
    % Get information from the dimensions of the input data
    num_worms = dims(1);
    final_t = dims(4);
    
    for curr_cluster_width = 1:length(cluster_width_cutoffs)
        j = cluster_width_cutoffs(curr_cluster_width);
        
        % Generate a random list of frames to sample
        frames = randi(final_t, 1, num_samples);
        
        bigger_counts = zeros(num_worms, num_samples);
        z_store = zeros(num_worms - 1, num_samples);
        
        for frame=1:length(frames)
            t = frames(frame);
            
            % Access the data for the primary worm node, initially for the first frame
            data1 = data(:,1,:,t);
            
            % Initialise empty matrices to store location coordinates
            coords = zeros(num_worms,2);
            
            coords(:,1) = data1(:,:,1);
            coords(:,2) = data1(:,:,2);
            
            Y = pdist(coords);
            
            Z = linkage(Y, 'complete');
            
            T = cluster(Z,'cutoff',j, 'criterion', 'distance');
            num_clusters = max(T);
            
            % Calculating the average number of worms per cluster
            a = unique(T);
            out = [a,histc(T(:),a)];
            counts_for_clusters = out(:,2).*(out(:,2)>=pop_thresh);
            wpc_counts(frame) = mean(counts_for_clusters);
            
            perc_worms_in_clust(frame) = sum(counts_for_clusters)/num_worms;
            
            for q = 1:length(counts_for_clusters)
                bigger_counts(q, frame) = counts_for_clusters(q);
            end
            
            for p = 1:num_worms-1
                z_store(p,frame) = Z(p,3);
            end
            
        end
        
        lin_counts = reshape(bigger_counts, 1, (num_worms*num_samples));
        lin_counts = lin_counts(lin_counts ~= 0);
        
        lin_heights = reshape(z_store,1,((num_worms-1)*num_samples));
        
        cluster_edges = [0:1:num_worms];
        wpc_hist = histcounts(lin_counts, cluster_edges, 'Normalization', 'probability');
        
        branch_edges = [0:0.5:20];
        branch_hist = histcounts(lin_heights, branch_edges, ...
            'Normalization', 'probability');
        
    end
    
elseif format_in == 'experiment'
    % Analagous code for profucing same histcount outputs from the
    % experimental .hdf5 output format
    
    for curr_cluster_width = 1:length(cluster_width_cutoffs)
        j = cluster_width_cutoffs(curr_cluster_width);
        % Randomly draw the frames to sample
        frames = data{3};
        to_sample = randi([min(frames),max(frames)], 1, num_samples);
        
        bigger_counts = [];
        z_store = [];
        
        for i = 1:num_samples;
            f = to_sample(i);
            
            returned = find(frames==f);
            
            while length(returned) < 2;
                f = randi([min(frames),max(frames)],1);
                returned = find(frames==f);
            end
            
            num_worms = length(returned);
            coords = zeros(num_worms,2);
            
            coords(:,1) = data{1}(returned);
            coords(:,2) = data{2}(returned);
            
            % Make pixel xy coordinates informative by converting to mm
            pix2mm = 0.1/19.5;
            coords = coords.*pix2mm;
            
            Y = pdist(coords);
            
            Z = linkage(Y, 'complete');
            
            T = cluster(Z,'cutoff',j, 'criterion', 'distance');
            num_clusters = max(T);
            
            % Calculating the average number of worms per cluster
            a = unique(T);
            out = [a,histc(T(:),a)];
            counts_for_clusters = out(:,2).*(out(:,2)>=pop_thresh);
            wpc_counts(i) = mean(counts_for_clusters);
            
            perc_worms_in_clust(i) = sum(counts_for_clusters)/num_worms;
            
            for q = 1:length(counts_for_clusters)
                bigger_counts(end+1) = counts_for_clusters(q);
            end
            
            for p = 1:num_worms-1
                z_store(end+1) = Z(p,3);
            end
            
        end
        
        lin_counts = bigger_counts;
        lin_counts = lin_counts(lin_counts ~= 0);
        
        cluster_edges = [0:1:40];
        wpc_hist = histcounts(lin_counts, cluster_edges, 'Normalization', 'probability');
        
        branch_edges = [0:0.5:20];
        branch_hist = histcounts(z_store, branch_edges, ...
            'Normalization', 'probability');
        
    end
    
    
end
end
