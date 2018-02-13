%% Function to output mean spread of points and std of median position over time
function [sig_x, std_sig_x] ...
    = inf_positionalmoments(data, format, fraction_to_sample)

% Specify fraction of frames to sample
if nargin<3
    fraction_to_sample = 0.1; % specify the proportion of frames to be sampled
end

if strcmp(format,'simulation') || strcmp(format,'complexsim')||strcmp(format,'simulation-test')
    burn_in = 0.5; % specifies how much to ignore at the start of the simulation
    
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
    num_worms = dims(1);
    final_frame = dims(4);
    
    % Sample fraction of the frames in the video
    num_samples = round(final_frame * (1 - burn_in) * fraction_to_sample);
    sampled_frames = randi([round(burn_in*final_frame) final_frame],1,num_samples);
    
    % initialise matrix to store spread of points
    std_pos = zeros(num_samples,1);
    
    for frameCtr=1:num_samples
        thisFrame = sampled_frames(frameCtr);
        
        % Access the data for the tracked worm node(s)
        thisFrameData = data(:,round(mean(trackedNodes)),:,thisFrame);
        
        % Initialise empty matrices to store location coordinates
        coords = zeros(num_worms,2);
        
        coords(:,1) = thisFrameData(:,:,1);
        coords(:,2) = thisFrameData(:,:,2);
        
        std_pos(frameCtr) = sqrt(sum(var(coords)));
    end
    
elseif format == 'experiment'
    % Analagous code for producing same histcount outputs from the
    % experimental .hdf5 output format
    
    % Make pixel xy coordinates informative by converting to mm
    pix2mm = 0.1/19.5;
    
    % Randomly draw the frames to sample
    frames = data{3};
    num_samples = floor(length(unique(frames)) * fraction_to_sample);
    frames_sampled = randi([min(frames),max(frames)], 1, num_samples);
    
    % initialise matrix to store spread of points
    std_pos = zeros(num_samples,1);
    
    for frameCtr = 1:num_samples
        thisFrame = frames_sampled(frameCtr);
        
        thisFrame_logInd = find(frames==thisFrame);
        
        while length(thisFrame_logInd) < 2 % resample if less than two worms in frame
            thisFrame = randi([min(frames),max(frames)],1);
            frames_sampled(frameCtr) = thisFrame;
            thisFrame_logInd = find(frames==thisFrame);
        end
        
        num_worms = length(thisFrame_logInd);
        coords = zeros(num_worms,2);
        
        coords(:,1) = data{1}(thisFrame_logInd).*pix2mm;
        coords(:,2) = data{2}(thisFrame_logInd).*pix2mm;
 
        std_pos(frameCtr) = sqrt(sum(var(coords)));            
    end
   
end

sig_x = mean(std_pos);
std_sig_x = std(std_pos);

end
