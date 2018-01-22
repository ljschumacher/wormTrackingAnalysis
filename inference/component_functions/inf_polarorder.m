% Function to calculate the polar order parameter over many sampled frames
function polar_order = inf_polarorder(data, format, fraction_to_sample)

if nargin<3
    fraction_to_sample = 0.1; % specify the proportion of frames to be sampled
end

if strcmp(format,'simulation') || strcmp(format,'complexsim')||strcmp(format,'simulation-test')
    burn_in = 0.25; % specifies how much to ignore at the start of the simulation
    L = 7.5;
    % Get the dimensions of the dataframe
    dims = size(data);
    if strcmp(format,'simulation')
        trackedNodes = 1:3;% only track nodes equivalent to the head
    elseif strcmp(format,'complexsim')
        trackedNodes = 1:6;
    elseif strcmp(format,'simulation-test')
        trackedNodes = 1;
    end
    
    % Get information from the dimensions of the input data
    final_frame = dims(4);
    % Sample fraction of the frames in the video
    num_samples = round(final_frame * (1 - burn_in) * fraction_to_sample);
    sampled_frames = randi([round(burn_in*final_frame) final_frame],1,num_samples); % could also sample without replacement, or regularly, but it doesn't seem to make a difference
    p_order_store = zeros(num_samples,1);
    % calculate orientations
    dxds = diff(data(:,fliplr(trackedNodes),1,:),1,2);
    dyds = diff(data(:,fliplr(trackedNodes),2,:),1,2);
    % correct for periodic boundary conditions
    dxds = correctForPeriodicBoundary(dxds,L);
    dyds = correctForPeriodicBoundary(dyds,L);
    % convert to orientation vectors
    orientations_x = squeeze(mean(dxds,2));
    orientations_y = squeeze(mean(dyds,2));
    % calculate orientation angles
    phis = atan2(orientations_y,orientations_x);
    for sampleCtr = 1:num_samples
        % calculate polar order parameter
        p_order_store(sampleCtr) = abs(mean(exp(1i*phis(:,sampled_frames(sampleCtr)))));
    end
    
    % Compute the average g(r) over the sampled timepoints
    polar_order = mean(p_order_store);
    
elseif format == 'experiment'
    % Analagous code for obtaining the same gr output from the
    % experimental .hdf5 data structure
    
    % Make pixel xy coordinates informative by converting to mm
    pix2mm = 0.1/19.5;
    
    % Randomly sample fraction of the frames in the video
    frames = data{3};
    num_samples = floor(length(unique(frames)) * fraction_to_sample);
    frames_sampled = randi([min(frames),max(frames)], 1, num_samples);
    
    % calculate orientations
    skelData = data{4};
    orientations_x = double(squeeze(skelData(1,1,:) - skelData(1,2,:)));
    orientations_y = double(squeeze(skelData(2,1,:) - skelData(2,2,:)));
    phis = atan2(orientations_y,orientations_x);
    p_order_store = zeros(num_samples,1);
    for sampleCtr = 1:num_samples
        thisFrame = frames_sampled(sampleCtr);
        
        thisFrame_idx = find(frames==thisFrame);
        
        while length(thisFrame_idx) < 2 % resample if less than two worms in frame
            thisFrame = randi([min(frames),max(frames)],1);
            frames_sampled(sampleCtr) = thisFrame;
            thisFrame_idx = find(frames==thisFrame);
        end
        p_order_store(sampleCtr) = abs(nanmean(exp(1i*phis(thisFrame_idx)))); % nanmean is necessary here as some objects are not skeletonised
    end
    
    % Compute the average over the sampled timepoints
    polar_order = mean(p_order_store);
end
end

function corrected = correctForPeriodicBoundary(raw, L)
corrected = raw;
corrected(raw<=-L/2) = raw(raw<=-L/2) + L;
corrected(raw>=L/2) = raw(raw>=L/2) - L;
end