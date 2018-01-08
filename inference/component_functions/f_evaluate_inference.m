function predictions_to_real = f_evaluate_inference(chosen_samples,...
    chosen_params,blindsim_dists,blind_file_names,hidden_params,test_params,...
    showPlot)

%function to evaluate the quality of inference by comparing predictions
%with hidden test samples

% For storing the predicted and real parameter values for the blind simulations
parameter_store = zeros(length(blind_file_names), length(test_params),2);

% For storing the log ratio of real to predicted parameters, along side the
% confidence interval the real paramater falls into, based on the posterioir
% distribution of the prediction
predictions_to_real = zeros(length(blind_file_names), length(test_params),2);

for blind_trial_ctr = 1:length(blind_file_names)
    %% use inference-returned simulation samples to make predictions for blind parameter values
    % construct weights for averaging samples
    chosen_samples_thistrial = chosen_samples(blind_trial_ctr,:,end); % use only samples for the smalles cut-off
    weights = 1./blindsim_dists(blind_trial_ctr,chosen_samples_thistrial,1);
    
%     % Report the closest single simulation match
%     single_best_match = chosen_samples_thistrial(1);
%     single_best_params = chosen_params(blind_trial_ctr,1,:);
    
    weighted_predict = sum(chosen_params(blind_trial_ctr,1,:).*weights)./sum(weights);
    
    % For each parameter being fitted, get the log ratio and CI
    for paramCtr = 1:length(test_params)
        predicted = weighted_predict(paramCtr);
        real_value = eval(['hidden_params.paramSamples.' test_params{paramCtr} '(blind_trial_ctr)']);
        difference = abs(predicted-real_value);
        
        % Store the raw parameter values
        parameter_store(blind_trial_ctr, paramCtr,1) = predicted;
        parameter_store(blind_trial_ctr, paramCtr,2) = real_value;
        
        %Compare the real value to the predicted value with a log ratio
        predictions_to_real(blind_trial_ctr, paramCtr,1) = ...
            abs(log(predicted/real_value));
        % Also report the associated confidence score, the interval around
        % the predicted value in which the real value lies
        this_chosen_params = chosen_params(blind_trial_ctr,:,paramCtr);
        interval_param_samples = this_chosen_params(this_chosen_params > (predicted-difference) ...
            & this_chosen_params < (predicted+difference));
        
        predictions_to_real(blind_trial_ctr, paramCtr,2) = ...
            1 - numel(interval_param_samples)/numel(this_chosen_params);
    end
end

if showPlot
    figure;
    for paramCtr = 1:length(test_params)
        subplot(1,length(test_params),paramCtr)
        scatter(parameter_store(:,paramCtr,1), parameter_store(:,paramCtr,2))
        hold on
        refline(1,0)
        xlabel('Predicted value')
        ylabel('Real value')
        title(test_params{paramCtr})
    end
end
end

