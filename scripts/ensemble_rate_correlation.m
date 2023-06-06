function [ensemble_rate_corr_across_sess,ensemble_rate_corr_elapsed_time] = ensemble_rate_correlation(activity_vectors_per_sess,multi_maps_mat,flag)

num_sessions = size(activity_vectors_per_sess,2);

if strcmp(flag,'full')
    ensemble_rate_vector_per_sess = [];
    for sess_ind = 1:num_sessions
        current_sess_activity = activity_vectors_per_sess{sess_ind};
        current_sess_ensemble_rate_vector = mean(current_sess_activity,'omitnan');
        
        ensemble_rate_vector_per_sess(:,sess_ind) = current_sess_ensemble_rate_vector;
    end
    
elseif strcmp(flag,'corr')
    ensemble_rate_vector_per_sess = activity_vectors_per_sess;
end

% calculating the ensemble rate correlation between pairs of sessions
ensemble_rate_corr_raw = corr(ensemble_rate_vector_per_sess,'rows','pairwise');

% Correcting for multiple maps - including data from pairs of
% sessions with the same map:
ensemble_rate_corr_across_sess = ensemble_rate_corr_raw .* multi_maps_mat;

% calculating ensemble rate correlation as function of elapsed time
ensemble_rate_corr_elapsed_time = [];
for time_interval = 1:num_sessions-1
    current_interval_values = diag(ensemble_rate_corr_across_sess,time_interval);
    ensemble_rate_corr_elapsed_time(time_interval) = mean(current_interval_values,'omitnan');
end

end
