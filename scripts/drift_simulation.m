function [full_drift_model,tuning_drift_model,rate_drift_model] = drift_simulation(simulation_parameters,state)

% loading random seed
rng(state)

% retrieving parameters
number_of_cells = simulation_parameters.number_of_cells;
tuning_smoothing_sigma = simulation_parameters.tuning_smoothing_sigma;
tuning_change_sigma = simulation_parameters.tuning_change_sigma;
activity_rate_mu = simulation_parameters.activity_rate_mu;
activity_rate_sigma = simulation_parameters.activity_rate_sigma;
activity_sigma_factor = simulation_parameters.activity_sigma_factor;
prob_tuning_change = simulation_parameters.prob_tuning_change;
prob_rate_change = simulation_parameters.prob_rate_change;
fraction_of_place_cells = simulation_parameters.fraction_of_place_cells;
number_of_bins = simulation_parameters.number_of_bins;
number_of_days = simulation_parameters.number_of_days;
number_of_simulations = simulation_parameters.number_of_simulations;
multi_maps_mat = ones(number_of_days);

full_drift_model = [];
tuning_drift_model = [];
rate_drift_model = [];

% starting a new loading bar
for simulation_ind = 1:number_of_simulations
     clc; 
     disp('Simulating drifting neuronal responses:')
     disp(['Realization: ',num2str(simulation_ind),'\',num2str(number_of_simulations)])
  
    % sample the activity rate of each neuron on
    % day 1 from a random log-normal distribution
    pd = makedist('Lognormal','mu',activity_rate_mu ,'sigma',activity_rate_sigma);
    activity_rate_for_each_neuron_day1 = random(pd,number_of_cells,1);
    
    % define a random tuning curve for each cell on day 1
    population_peak_locations_day1 = [];
    for cell_ind = 1:number_of_cells
        current_cell_peak_location = randperm(number_of_bins,1);
        
        current_cell_tuning_curve = zeros(1,number_of_bins);
        current_cell_tuning_curve(current_cell_peak_location) = 1;
        
        current_cell_tuning_curve_smoothed = imgaussfilt(current_cell_tuning_curve,tuning_smoothing_sigma).*activity_rate_for_each_neuron_day1(cell_ind);
        
        population_peak_locations_day1(cell_ind,1) = current_cell_peak_location;
    end
    
    
    % calculating the change in activity and tuning of each cell across days
    positional_shifts = [];
    activity_shifts = [];
    for day_ind = 2:number_of_days
        % calculate the change in tuning for each day
        tuning_change_binary = rand(number_of_cells,1) <= prob_tuning_change;
        pd = makedist('Normal','mu',0,'sigma',tuning_change_sigma);
        t = truncate(pd,-10,10);
        positional_shifts(:,day_ind) = round(random(t,number_of_cells,1)) .* tuning_change_binary;
        
        % calculate the change in activity for each day
        for cell_ind = 1:number_of_cells
            current_cell_mu = activity_rate_for_each_neuron_day1(cell_ind);
            current_cell_sigma = current_cell_mu./activity_sigma_factor;
            pd = makedist('Normal','mu',0,'sigma',current_cell_sigma);
            t = truncate(pd,-current_cell_sigma*3,current_cell_sigma*3);
            activity_shifts(cell_ind,day_ind) = random(t,1,1);
        end
    end
    
    
    % calculating the drift in representations over days
    population_vectors_across_days = [];
    population_vectors_across_days_no_tuning_shift = [];
    population_vectors_across_days_no_rate_shift = [];
    
    ensemble_rate_vectors_across_days = [];
    ensemble_rate_vectors_across_days_no_tuning_shift = [];
    ensemble_rate_vectors_across_days_no_rate_shift = [];
    
    
    PC_ind_each_day = {};
    PC_ind_each_day_no_rate_shift = {};
    for day_ind = 1:number_of_days
        
        PCs_each_day = ones(number_of_cells,1);
        
        if day_ind == 1 % day 1 does not change
            current_day_activity = activity_rate_for_each_neuron_day1;
            current_day_peak_locations = population_peak_locations_day1;
        else
            current_day_activity = current_day_activity_yesterday + activity_shifts(:,day_ind);
            current_day_peak_locations = current_day_peak_locations_yesterday + positional_shifts(:,day_ind);
            
            % correcting positional shifts
            above_boundary = find(current_day_peak_locations>number_of_bins);
            corrected_peak_location_above = (number_of_bins+1)-(current_day_peak_locations(above_boundary)-number_of_bins);
            
            current_day_peak_locations(above_boundary) = corrected_peak_location_above;
            
            below_boundary = find(current_day_peak_locations<1);
            corrected_peak_location_below = (1-current_day_peak_locations(below_boundary));
            current_day_peak_locations(below_boundary) = corrected_peak_location_below;
        end
        
        population_tuning_curves_smoothed = [];
        population_tuning_curves_smoothed_no_tuning_shift = [];
        population_tuning_curves_smoothed_no_rate_shift = [];
        for cell_ind = 1:number_of_cells
            current_cell_peak_location = current_day_peak_locations(cell_ind); % full and tuning drift model
            current_cell_peak_location_no_shift = population_peak_locations_day1(cell_ind); % rate drift model
            
            if rand(1) <= fraction_of_place_cells % if current cell is  not aplace cell
                PCs_each_day(cell_ind) = 0;
                
                current_cell_tuning_curve = rand(1,number_of_bins).*0.9; %0.8 % defining a new noisy tuning curve
                
                % deining peak location for rate drift model
                current_cell_tuning_curve_no_shift = current_cell_tuning_curve;
                current_cell_tuning_curve_no_shift(current_cell_peak_location_no_shift) = 1;
                current_cell_tuning_curve_no_shift = current_cell_tuning_curve_no_shift./sum(current_cell_tuning_curve_no_shift);
                
                
                % deining peak location for full and tuning drift model
                current_cell_tuning_curve(current_cell_peak_location) = 1;
                current_cell_tuning_curve = current_cell_tuning_curve./sum(current_cell_tuning_curve);
                
            else
                current_cell_tuning_curve = zeros(1,number_of_bins);
                current_cell_tuning_curve(current_cell_peak_location) = 1;
                current_cell_tuning_curve = current_cell_tuning_curve./sum(current_cell_tuning_curve);
                
                % deining peak location for rate drift model
                current_cell_tuning_curve_no_shift = zeros(1,number_of_bins);
                current_cell_tuning_curve_no_shift(current_cell_peak_location_no_shift) = 1;
                current_cell_tuning_curve_no_shift = current_cell_tuning_curve_no_shift./sum(current_cell_tuning_curve_no_shift);
            end
            
            
            current_cell_tuning_curve_smoothed = imgaussfilt(current_cell_tuning_curve,tuning_smoothing_sigma);
            current_cell_tuning_curve_smoothed_no_shift = imgaussfilt(current_cell_tuning_curve_no_shift,tuning_smoothing_sigma);
            
            % store smoothed tuning curves for all cells
            population_tuning_curves_smoothed(cell_ind,:) = current_cell_tuning_curve_smoothed.*current_day_activity(cell_ind);
            population_tuning_curves_smoothed_no_tuning_shift(cell_ind,:) = current_cell_tuning_curve_smoothed_no_shift.*current_day_activity(cell_ind);
            population_tuning_curves_smoothed_no_rate_shift(cell_ind,:) = current_cell_tuning_curve_smoothed.*activity_rate_for_each_neuron_day1(cell_ind);
        end
        
        PC_ind_each_day_no_rate_shift{day_ind} = find(PCs_each_day)';
        
        % correcting activity shifts that lead to negative activity rates
        below_rate_boundary = current_day_activity > 0;
        rate_change_binary = rand(number_of_cells,1) <= prob_rate_change;
        
        peak_locations_corrected = current_day_peak_locations.*(rate_change_binary & below_rate_boundary & PCs_each_day);
        peak_locations_corrected(peak_locations_corrected==0) = NaN;
        
        % if a cell became not active than it is no longer a place cell
        PCs_each_day(~(rate_change_binary & below_rate_boundary)) = 0;
        PC_ind_each_day{day_ind} = find(PCs_each_day)';
        
        
        population_vectors_across_days(:,:,day_ind) = population_tuning_curves_smoothed.*(rate_change_binary& below_rate_boundary);
        population_vectors_across_days_no_tuning_shift(:,:,day_ind) = population_tuning_curves_smoothed_no_tuning_shift.*(rate_change_binary& below_rate_boundary);
        population_vectors_across_days_no_rate_shift(:,:,day_ind) = population_tuning_curves_smoothed_no_rate_shift.*(rate_change_binary);
        
        ensemble_rate_vectors_across_days(:,day_ind) = current_day_activity.*(rate_change_binary & below_rate_boundary);
        ensemble_rate_vectors_across_days_no_tuning_shift(:,day_ind) = current_day_activity.*(rate_change_binary & below_rate_boundary);
        ensemble_rate_vectors_across_days_no_rate_shift(:,day_ind) = activity_rate_for_each_neuron_day1.*(rate_change_binary);
        
        current_day_peak_locations_yesterday = current_day_peak_locations;
        current_day_activity_yesterday = current_day_activity;
        
    end
    
    
    % Calculate the across days stability measurements for the "Full drift" model:
    % Population vector correlation
    [pv_corr_mat_across_days,mean_pv_corr_across_days,mean_pv_corr_elapsed_time] = population_vector_correlation(population_vectors_across_days,PC_ind_each_day,multi_maps_mat);
    full_drift_model.population_vectors_each_day{simulation_ind} = population_vectors_across_days;
    full_drift_model.PC_ind_each_day(simulation_ind,:) = PC_ind_each_day;
    full_drift_model.pv_corr_mat(:,:,simulation_ind) = pv_corr_mat_across_days;
    %   full_drift_model.mean_pv_corr(:,:,simulation_ind) = mean_pv_corr_across_days;
    full_drift_model.pv_corr_elapsed_time(simulation_ind,:) = mean_pv_corr_elapsed_time;
    
    % Ensemble rate correlation
    [~,ensemble_rate_corr_elapsed_time] = ensemble_rate_correlation(ensemble_rate_vectors_across_days,multi_maps_mat,'corr');
    full_drift_model.ensemble_rate_corr_elapsed_time(simulation_ind,:) = ensemble_rate_corr_elapsed_time;
    
    % Tuning curve correlation
    [~,mean_tuning_curve_corr_elapsed_time] = tuning_curve_correlation(population_vectors_across_days,PC_ind_each_day,multi_maps_mat);
    full_drift_model.tuning_curve_corr_elapsed_time(simulation_ind,:) = mean_tuning_curve_corr_elapsed_time;
    
    % Population vector template matching decoder
    [mean_decoder_accuracy_elapsed_time] = population_vector_decoder(population_vectors_across_days,PC_ind_each_day,multi_maps_mat);
    full_drift_model.decoder_accuracy_elapsed_time(simulation_ind,:) = mean_decoder_accuracy_elapsed_time;
    
    
    % Calculate the across days stability measurements for the "Rate drift" model:
    % Population vector correlation
    [pv_corr_mat_across_days,mean_pv_corr_across_days,mean_pv_corr_elapsed_time] = population_vector_correlation(population_vectors_across_days_no_tuning_shift,PC_ind_each_day,multi_maps_mat);
    rate_drift_model.pv_corr_mat(:,:,simulation_ind) = pv_corr_mat_across_days;
    rate_drift_model.mean_pv_corr(:,:,simulation_ind) = mean_pv_corr_across_days;
    rate_drift_model.pv_corr_elapsed_time(simulation_ind,:) = mean_pv_corr_elapsed_time;
    
    % Ensemble rate correlation
    [~,ensemble_rate_corr_elapsed_time] = ensemble_rate_correlation(ensemble_rate_vectors_across_days_no_tuning_shift,multi_maps_mat,'corr');
    rate_drift_model.ensemble_rate_corr_elapsed_time(simulation_ind,:) = ensemble_rate_corr_elapsed_time;
    
    % Tuning curve correlation
    [~,mean_tuning_curve_corr_elapsed_time] = tuning_curve_correlation(population_vectors_across_days_no_tuning_shift,PC_ind_each_day,multi_maps_mat);
    rate_drift_model.tuning_curve_corr_elapsed_time(simulation_ind,:) = mean_tuning_curve_corr_elapsed_time;
    
    % Population vector template matching decoder
    [mean_decoder_accuracy_elapsed_time] = population_vector_decoder(population_vectors_across_days_no_tuning_shift,PC_ind_each_day,multi_maps_mat);
    rate_drift_model.decoder_accuracy_elapsed_time(simulation_ind,:) = mean_decoder_accuracy_elapsed_time;
    
    % Calculate the across days stability measurements for the "Tuning drift" model:
    % Population vector correlation
    [pv_corr_mat_across_days,mean_pv_corr_across_days,mean_pv_corr_elapsed_time] = population_vector_correlation(population_vectors_across_days_no_rate_shift,PC_ind_each_day,multi_maps_mat);
    tuning_drift_model.pv_corr_mat(:,:,simulation_ind) = pv_corr_mat_across_days;
    tuning_drift_model.mean_pv_corr(:,:,simulation_ind) = mean_pv_corr_across_days;
    tuning_drift_model.pv_corr_elapsed_time(simulation_ind,:) = mean_pv_corr_elapsed_time;
    
    % Ensemble rate correlation
    [~,ensemble_rate_corr_elapsed_time] = ensemble_rate_correlation(ensemble_rate_vectors_across_days_no_rate_shift,multi_maps_mat,'corr');
    tuning_drift_model.ensemble_rate_corr_elapsed_time(simulation_ind,:) = ensemble_rate_corr_elapsed_time;
    
    % Tuning curve correlation
    [~,mean_tuning_curve_corr_elapsed_time] = tuning_curve_correlation(population_vectors_across_days_no_rate_shift,PC_ind_each_day,multi_maps_mat);
    tuning_drift_model.tuning_curve_corr_elapsed_time(simulation_ind,:) = mean_tuning_curve_corr_elapsed_time;
    
    % Population vector template matching decoder
    [mean_decoder_accuracy_elapsed_time] = population_vector_decoder(population_vectors_across_days_no_rate_shift,PC_ind_each_day,multi_maps_mat);
    tuning_drift_model.decoder_accuracy_elapsed_time(simulation_ind,:) = mean_decoder_accuracy_elapsed_time;
    
end

end

