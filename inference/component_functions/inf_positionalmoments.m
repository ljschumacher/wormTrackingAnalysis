%% Function to output mean spread of points and std of median position over time
function [sig_x, kurt_x] ...
    = inf_positionalmoments(data, format, fraction_to_sample)

% Specify fraction of frames to sample
if nargin<3
    fraction_to_sample = 0.1; % specify the proportion of frames to be sampled
end

if strcmp(format,'simulation') || strcmp(format,'complexsim')||strcmp(format,'simulation-test')
    burn_in = 0.5; % specifies how much to ignore at the start of the simulation
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
    num_worms = dims(1);
    final_frame = dims(4);
    
    % Sample fraction of the frames in the video
    num_samples = round(final_frame * (1 - burn_in) * fraction_to_sample);
    sampled_frames = randi([round(burn_in*final_frame) final_frame],1,num_samples);
    
    % initialise matrix to store spread of points
    std_pos = zeros(num_samples,1);
    kurt_pos = zeros(num_samples,1);
    
    for frameCtr=1:num_samples
        thisFrame = sampled_frames(frameCtr);
        
        % Access the data for the tracked worm node(s)
        thisFrameData = data(:,round(mean(trackedNodes)),:,thisFrame);
        
        % Initialise empty matrices to store location coordinates
        coords = zeros(num_worms,2);
        
        coords(:,1) = thisFrameData(:,:,1);
        coords(:,2) = thisFrameData(:,:,2);
        
        % center each frame on center of mass - needed for periodic boundary conditions
        % this calculates the centre of mass for periodic boundaries
        c_x = mean(cos(coords(:,1)/L*2*pi));
        s_x = mean(sin(coords(:,1)/L*2*pi));
        xmean = L/2/pi*(atan2(-s_x,-c_x) + pi);
        c_y = mean(cos(coords(:,2)/L*2*pi));
        s_y = mean(sin(coords(:,2)/L*2*pi));
        ymean = L/2/pi*(atan2(-s_y,-c_y) + pi);
        
        xdistanceFromMean = periodiceucdist(coords(:,1),xmean);
        ydistanceFromMean = periodiceucdist(coords(:,2),ymean);
        thisN = size(coords,1);
        varPbc = [sum(xdistanceFromMean.^2), sum(ydistanceFromMean.^2)]./(thisN-1);
        std_pos(frameCtr) = sqrt(sum(varPbc));
        
        kurtx = mean(xdistanceFromMean.^4)./(mean(xdistanceFromMean.^2)).^2;
        kurty = mean(ydistanceFromMean.^4)./(mean(ydistanceFromMean.^2)).^2;
        kurt_pos(frameCtr) = mean([kurtx kurty]);
        
        %         c0 = coords - mean(coords);
        %         rad_gyr(frameCtr) = sum(sum(c0.*c0))/num_worms; % this gives the same as the variane above
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
    kurt_pos = zeros(num_samples,1);
    
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
        
        std_pos(frameCtr) = sqrt(sum(var(coords))); % this is like the radius of gyration
        kurt_pos(frameCtr) = mean(kurtosis(coords,0));
    end
    
end

sig_x = mean(std_pos);
% std_sig_x = std(std_pos);
kurt_x = mean(kurt_pos);
% std_kurt_x = std(kurt_pos);

end
