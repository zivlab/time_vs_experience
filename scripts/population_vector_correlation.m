function [pv_corr_mat_across_sess,mean_pv_corr_across_sess,mean_pv_corr_elapsed_time] = population_vector_correlation(population_vectors_per_sess,place_cells_indices_per_sess,multi_maps_mat)

num_sessions = size(multi_maps_mat,1);

pv_corr_mat_across_sess = {};
mean_pv_corr_across_sess = [];
for sess1_ind = 1:num_sessions
    population_vectors_sess1 = population_vectors_per_sess(:,:,sess1_ind);
    place_cells_ind_sess1 = place_cells_indices_per_sess{sess1_ind};
    
    for sess2_ind = 1:num_sessions
        population_vectors_sess2 = population_vectors_per_sess(:,:,sess2_ind);
        place_cells_ind_sess2 = place_cells_indices_per_sess{sess2_ind};
        
        % subset cells that were found to be place cells in at least one
        % of the two sessions:
        place_cells_atleast_one_sess = unique([place_cells_ind_sess1,place_cells_ind_sess2]);
        
        pv_only_place_cells_sess1 = population_vectors_sess1(place_cells_atleast_one_sess,:);
        pv_only_place_cells_sess2 = population_vectors_sess2(place_cells_atleast_one_sess,:);
        
        % subset cells that were found to be active in both sessions (at
        % least 1 calcium event):
        active_cells_ind_sess1 = mean( pv_only_place_cells_sess1,2,'omitnan') > 0;
        active_cells_ind_sess2 = mean( pv_only_place_cells_sess2,2,'omitnan') > 0;
        active_cells_ind_both_sess = active_cells_ind_sess1 & active_cells_ind_sess2;
        
        pv_only_active_place_cells_sess1 = pv_only_place_cells_sess1(active_cells_ind_both_sess,:);
        pv_only_active_place_cells_sess2 = pv_only_place_cells_sess2(active_cells_ind_both_sess,:);
        
         % Correcting for multiple maps - including data from pairs of
        % sessions with the same map:
        if ~isnan(multi_maps_mat(sess1_ind,sess2_ind)) && sum(active_cells_ind_both_sess) > 0 % same map and check overlap in active cells
            
            % calculating population vector correlation
            pv_corr_mat = corr(pv_only_active_place_cells_sess1,pv_only_active_place_cells_sess2,'rows','pairwise');
            
            pv_corr_mat_across_sess{sess1_ind,sess2_ind} = pv_corr_mat;
            mean_pv_corr_across_sess(sess1_ind,sess2_ind) = mean(diag(pv_corr_mat),'omitnan');
            
        else % different maps
            pv_corr_mat_across_sess{sess1_ind,sess2_ind} = nan(size(pv_corr_mat,1));
            mean_pv_corr_across_sess(sess1_ind,sess2_ind) = NaN;
        end
    end
end

% combine all pv corr matrices to create a single pv corr matrix across days
pv_corr_mat_across_sess = cell2mat(pv_corr_mat_across_sess);

% calculating pv correlation as function of elapsed time
mean_pv_corr_elapsed_time = [];
for time_interval = 1:num_sessions-1
    current_interval_values = diag(mean_pv_corr_across_sess,time_interval);
    mean_pv_corr_elapsed_time(time_interval) = mean(current_interval_values,'omitnan');
end


end



