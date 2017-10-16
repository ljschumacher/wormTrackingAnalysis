% Function to calculate the pair correlation function over many sampled
% frames, and return a discretised array of the resulting g(r) distribution
function gr3 = inf_gr(data, format)

% Create bins of a given width, to store the data in. The bin width
% controls the width of the rings drawn from each reference particle during
% the computation of g(r). It also controls the bins into which the g(r)
% distribution is discretised.
bin_width = 0.125;
bins = 0:bin_width:10;

%Also specify the proportion of frames to be sampled e.g to sample 20% of
%the frames, use the following: 'to_sample = 0.2'.
to_sample = 0.2;

if format == 'simulation'
    
    % Get the dimensions of the dataframe
    dims = size(data);
    data = data(:,1,:,:);
    
    % Get information from the dimensions of the input data
    num_worms = dims(1);
    final_t = dims(4);
    
    % Sample 20% of the frames in the video
    num_samples = floor(final_t * to_sample);
    sampled_t = randi([1 final_t],1,num_samples);
    
    for t = 1:num_samples
        
        % Access the data for the primary worm node, initially for the first frame
        data1 = data(:,1,:,sampled_t(t));
        
        % Initialise empty matrices to store location coordinates
        coords = zeros(num_worms,2);
        
        coords(:,1) = data1(:,:,1);
        coords(:,2) = data1(:,:,2);
        
        % Calculate pairwise distances with custom distance function
        % 'periodiceucdists' to take into account the horizontal and
        % vertical periodicity of the simulations.
        pair_dist = pdist(coords, @periodiceucdist);
        
        % Get the histogram counts of the pair_dist data using the bins
        gr1 = histcounts(pair_dist,bins,'Normalization','probability');
        
        % Radial distribution function
        % Normalization step
        R = max(bins);
        gr2 = gr1.*R^2./(2*bins(2:end)*bin_width*((num_worms^2)-num_worms));
        
        % Store the gr information for each of the sampled timepoints
        if t == 1
            gr_store = zeros(num_samples,length(gr2));
        end
        gr_store(t,:) = gr2;
    end
    
    % Compute the average g(r) over the sampled timepoints
    gr3 = zeros(1,length(gr2));
    
    for i = 1:length(gr2)
        gr3(i) = mean(gr_store(:,i));
    end
    
elseif format == 'complexsim'
    
    % Get the dimensions of the dataframe
    dims = size(data);
    data = data(:,:,:,:);
    
    % Get information from the dimensions of the input data
    num_worms = dims(1);
    final_t = dims(4);
    
    % Sample 20% of the frames in the video
    num_samples = floor(final_t * 0.2);
    sampled_t = randi([1 final_t],1,num_samples);
    
    for t = 1:num_samples
        
        % Access the data for the 3rd worm node, initially for the first frame
        data1 = data(:,3,:,sampled_t(t));
        
        % Initialise empty matrices to store location coordinates
        coords = zeros(num_worms,2);
        
        coords(:,1) = data1(:,:,1);
        coords(:,2) = data1(:,:,2);
        
        % Calculate pairwise distances with custom distance function
        % 'periodiceucdists' to take into account the horizontal and
        % vertical periodicity of the simulations.
        pair_dist = pdist(coords, @periodiceucdist);
        
        % Get the histogram counts of the pair_dist data using the bins
        gr1 = histcounts(pair_dist,bins,'Normalization','probability');
        
        % Radial distribution function
        % Normalization step
        R = max(bins);
        gr2 = gr1.*R^2./(2*bins(2:end)*bin_width*((num_worms^2)-num_worms));
        
        % Store the gr information for each of the sampled timepoints
        if t == 1
            gr_store = zeros(num_samples,length(gr2));
        end
        gr_store(t,:) = gr2;
    end
    
    % Compute the average g(r) over the sampled timepoints
    gr3 = zeros(1,length(gr2));
    
    for i = 1:length(gr2)
        gr3(i) = mean(gr_store(:,i));
    end
    
elseif format == 'experiment'
    % Analagous code for obtaining the same gr output from the
    % experimental .hdf5 data structure
    frames = data{3};
    
    % Randomly sample 20% of the frames in the video
    num_samples = floor(length(unique(frames)) * 0.2);
    to_sample = randi([min(frames),max(frames)], 1, num_samples);
    
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
        
        % Obtain the pairwise distances with pdist
        pair_dist = pdist(coords);
        
        % Get the histogram counts of the pair_dist data using the bins
        gr1 = histcounts(pair_dist,bins,'Normalization','probability');
        
        % Radial distribution function
        % Normalization step
        R = max(bins);
        gr2 = gr1.*R^2./(2*bins(2:end)*bin_width*((num_worms^2)-num_worms));
        
        % Store the gr information for each of the sampled timepoints
        if i == 1
            gr_store = zeros(num_samples,length(gr2));
        end
        gr_store(i,:) = gr2;
    end
    
    % Compute the average g(r) over the sampled timepoints
    gr3 = zeros(1,length(gr2));
    
    for p = 1:length(gr2)
        gr3(p) = mean(gr_store(:,p));
    end   
end
end