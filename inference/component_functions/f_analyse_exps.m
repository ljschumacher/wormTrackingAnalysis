function [exp_ss_array, strain_list] = f_analyse_exps(strain_list, dataset,...
    num_statistics)

addpath('../')

useJoinedTraj = true;

% Construct arrays for storing the summary statistic outputs
exp_ss_array = cell(length(strain_list), 1 + num_statistics);

% set the filtering parameters
if dataset == 1
    intensityThreshold = 50;
elseif dataset == 2
    intensityThreshold = 60;
end
maxBlobSize = 1e4; % chosen large as to not filter out blobs of multiple pharynxes, as we may still want to count these in pair correlation and hierarchical clustering
pixelsize = 100/19.5; % 100 microns are 19.5 pixels;

% For each one of these list files provided i.e. for each group of videos..
for strainCtr = 1:length(strain_list)
    %% load file lists
    if dataset == 1
        [phaseFrames,filenames,~] = xlsread(['../datalists/' strain_list{strainCtr} '_40_list.xlsx'],1,'A1:E15','basic');
    elseif dataset == 2
        [phaseFrames,filenames,~] = xlsread(['../datalists/' strain_list{strainCtr} '_40_g_list.xlsx'],1,'A1:E15','basic');
    end
    if ~useJoinedTraj
        filenames = strrep(filenames,'/data2/shared/data/twoColour/Results/',...
            '/end/home/lschumac/databackup/data/twoColour/ResultsUnjoinedTrajectories/');
    end
    
    num_expmnts = length(filenames);
    
    exp_replicate_ss_array = cell(num_expmnts, 1 + num_statistics);
    %% For each movie file identified in this list...
    for expCtr = 1:num_expmnts
        filename = filenames{expCtr};
        
        %% OBTAINING AND FILTERING DATA
        
        % Must obtain logical indices for filtering the tracking data based
        % on the intensity and size of potential worms. This is to
        % eliminate larvae from the files
        
        % Read in the trajectory data
        trajData = h5read(filename, '/trajectories_data');
        % Read in the associated data
        blobData = h5read(filename, '/blob_features');
        
        % Get filtering indices according to the parameters outlined
        filter_logInd = filterIntensityAndSize(blobData,pixelsize,...
            intensityThreshold,maxBlobSize);%&has_skeleton;
        % filter for experimental phase
        [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,'joining',expCtr);
        phaseFilter_logInd = trajData.frame_number < lastFrame & trajData.frame_number > firstFrame;
        filter_logInd(~phaseFilter_logInd)=false;
        
        % Then split it to obtain the important information, whilst
        % applying the filter obtained from \blob_features above
        x_data = trajData.coord_x(filter_logInd);
        y_data = trajData.coord_y(filter_logInd);
        frames = trajData.frame_number(filter_logInd);
        %         relative_blob_sizes = blobData.area(filter_logInd)./mode(blobData.area(filter_logInd));
        if num_statistics>=4
            % also pass on skeleton data to compute orientational order
            skelData = h5read(filename,'/skeleton');
%             has_skeleton = squeeze(~any(any(isnan(skelData))));
            skelData = skelData(:,:,filter_logInd);
            in_data = {x_data, y_data, frames, skelData};
        else
            in_data = {x_data, y_data, frames};
        end
        
        % Then compute all chosen summary statistics, as with the simulated
        % data above.
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        fraction_to_sample = 1/max(frameRate,1)/3; % specifiy fraction of frames to sample
        ss_results = f_compute_ss(in_data, 'experiment', fraction_to_sample, num_statistics);
        
        for each_ss = 1:length(ss_results)
            exp_replicate_ss_array{expCtr, each_ss+1} = ss_results{each_ss};
        end
        % display progress
        disp(['analysed ' num2str(expCtr) '/' num2str(num_expmnts) ' experiments for strain ' strain_list{strainCtr}])
    end
    
    % Combine SS from different experimental replicates, to give single
    % reference for each strain
    replicate_summary = cell(1,num_statistics);
    for statCtr = 1:num_statistics
        replicate_summary{statCtr} = mean(cell2mat(exp_replicate_ss_array(:,1+statCtr)));
    end
    exp_ss_array{strainCtr, 1} = strain_list{strainCtr};
    exp_ss_array(strainCtr, 2:end) = replicate_summary;
end
end