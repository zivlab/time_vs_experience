%%  Defining basic paramaters and loading required data
clear; clc;

% Set paths
repository_path = 'Z:\personal\daniel\phd\time_vs_experience\neuron_published_dataset';

data_path = [repository_path,'\data\'];
scripts_path = [repository_path,'\scripts\'];
addpath(scripts_path)  

% Setting parameters
mouse_list = {'C82M2','C83M3','C84M2','C82M3','C82M4','C84M4','C89M0','C92M6'};
number_of_mice = length(mouse_list);
day_ind_per_env = {2:2:20,1:4:21};
elapsed_days_per_env = {2:2:18;4:4:20};
env_colors = [10 105 170; 190 25 50]./255;
env_full_names = {'Environment A','Environment B'};

% Load processed data and colormaps
load([data_path,'processed_data.mat'])
load([data_path,'example_fov.mat'])
load([data_path,'simulation_rng.mat'])
load([data_path,'magma_colormap.mat'])

%% Figure 2C - Cell registration demonstration for an example FOV
sess_list = [1,2,10];
day_list = [2,4,20];

figure('unit','normalized','position',[0.25 0.15 0.3 0.55])
for sess_ind = 1:length(sess_list)
    current_sess_fov = zeros(size(example_fov,1),size(example_fov,2),length(sess_list));
    current_sess_fov(:,:,sess_ind) = example_fov(:,:,sess_list(sess_ind));
    
    subplot(2,2,sess_ind)
    imagesc(current_sess_fov)
    set(gca,'xtick',[],'ytick',[])
    title(['Session ',num2str(sess_list(sess_ind)),' (day ',num2str(day_list(sess_ind)),'):'],'FontSize',14,'Fontweight','normal');
    axis square
end

fov_overlay = example_fov(:,:,sess_list);
subplot(2,2,4)
imagesc(fov_overlay)
set(gca,'xtick',[],'ytick',[])
title(['Overlay:'],'FontSize',14,'Fontweight','normal');
axis square

%% Figure 2D - Reponsiveness of all place cells across sessions

% IMPORTANT NOTE%
% Please note that it is required to save the plots as PDF of highest resolution
% in order to obtain the same visuals as in the published manuscript.

tuning_curves = population_vectors_not_smoothed_all_mice;
place_cells_ind = PC_ind_each_sess_both_sides_all_mice;
smoothing_sigma = 1.75;

% Visulizing neuronal responses for six sessions in each environment
zivplots(tuning_curves,place_cells_ind,smoothing_sigma)


%% Figure 2E - Population vector correlation across environments and running directions
mean_pv_corr_both_dirs_both_env_across_mice = mean(pv_corr_both_dirs_both_env_all_mice,3,'omitnan');

figure('units','normalized','position',[0.3 0.25 0.2 0.4])
imagesc(mean_pv_corr_both_dirs_both_env_across_mice,[0 0.55])
hold on
plot(xlim,[20.5 20.5],'k');plot([20.5 20.5],ylim,'k');plot(xlim,[60.5 60.5],'k');plot([60.5 60.5],ylim,'k');
plot(xlim,[40.5 40.5],'w','linewidth',1.5);plot([40.5 40.5],ylim,'w','linewidth',1.5)
title('Figure 2E','FontSize',24,'Fontweight','normal')
set(gca,'xtick',10:20:80,'xticklabels',{'Right','Left'},'ytick',10:20:80,'yticklabels',{'Right','Left'})
ytickangle(90)
colormap(jet)
axis square
cb = colorbar();
cb.Ticks = [0 0.55];
ylabel(cb,'PV correlation','FontSize',12,'Rotation',270);
ylabel('  Environment B       Environment A')
xlabel('Environment A      Environment B')

%% Figure 2F - Population vector correlation across sessions for each environment

figure('units','normalized','position',[0.15 0.2 0.675 0.5])
axes('pos',[0 1 1 1],'visible','off');
text(0.45,-0.05,'Figure 2F','fontsize',24);
for env_ind = 1:2
    current_env_days = day_ind_per_env{env_ind};
    current_env_pv_corr_mat = cell2mat(pv_corr_mat_all_mice(:,env_ind,:));
    avg_pv_corr_mat_across_mice = mean(current_env_pv_corr_mat,3,'omitnan');
    
    current_env_pv_corr_mat_inset = cell2mat(mean_pv_corr_all_mice(:,env_ind,:));
    avg_pv_corr_mat_across_mice_inset =  mean(current_env_pv_corr_mat_inset,3,'omitnan');
    avg_pv_corr_mat_across_mice_inset(logical(eye(size(avg_pv_corr_mat_across_mice_inset,1)))) = NaN;
    
    
    sub1 = axes();
    set(sub1,'units','normalized','position',[0.0715 0.15 0.275 0.65]+(env_ind-1)*[0.485 0 0 0])
    imagesc(avg_pv_corr_mat_across_mice,[-0.1 1])
    hold on
    for line = 1:(size(avg_pv_corr_mat_across_mice,1)./20)-1
        plot([20.5 20.5]+20*(line-1),ylim,'color','k')
        plot(xlim,[20.5 20.5]+20*(line-1),'color','k')
    end
    set(gca,'xtick',10.5:20:size(avg_pv_corr_mat_across_mice,1),'xticklabel',current_env_days,...
        'ytick',10.5:20:size(avg_pv_corr_mat_across_mice,1),'yticklabel',current_env_days)
    xlabel('Day','FontSize',16)
    ylabel('Day','FontSize',16)
    title([env_full_names{env_ind},':'],'units','normalized','Position',[0.22 1.01],'FontSize',16,'Fontweight','normal');
    cb = colorbar();
    cb.Position = [0.3575 0.15 0.015 0.325]+(env_ind-1)*[0.485 0 0 0];
    cb.Ticks = [-0.1 1];
    ylabel(cb,'PV correlation','FontSize',12,'Rotation',270);
    
    
    sub2 = axes();
    set(sub2,'units','normalized','position',[0.3575 0.55 0.11 0.25]+(env_ind-1)*[0.485 0 0 0])
    imagesc(avg_pv_corr_mat_across_mice_inset,'AlphaData',~isnan(avg_pv_corr_mat_across_mice_inset),[0 0.55])
    set(gca,'xtick',2:2:10,'xticklabel',current_env_days(2:2:end),'ytick',[],'color',0*[1 1 1])
    title({'Average across';'matching positions'},'units','normalized','Position',[0.5 1.01],'FontSize',12,'Fontweight','normal');
    cb = colorbar();
    cb.Position = [0.4775 0.55 0.0125 0.25]+(env_ind-1)*[0.485 0 0 0];
    cb.Ticks = [0 0.55];
    ylabel(cb,'PV correlation','FontSize',10,'Rotation',270);
    
    colormap(jet)
end

%% Figure 3B - Ensemble rate correlation as a function of elapsed time or number of sessions between visits
figure_name = 'Figure 3B';
stability_measurement = ensemble_rate_corr_elapsed_time_all_mice;
measurement_name = 'Ensemble rate correlation';
ylims = [0 0.7];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');

clc;
disp('Testing for the effects of time and experience on')
disp('representational drift in ensemble rate correlations')
disp('using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure 3C - Tuning curve correlation as a function of elapsed time or number of sessions between visits
figure_name = 'Figure 3C';
stability_measurement = tuning_curve_corr_elapsed_time_all_mice;
measurement_name = 'Tuning curve correlation';
ylims = [0 0.7];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');

clc;
disp('Testing for the effects of time and experience on')
disp('representational drift in tuning curve correlations')
disp('using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure 3D - Population vector correlation as a function of elapsed time or number of sessions between visits
figure_name = 'Figure 3D';
stability_measurement = pv_corr_elapsed_time_all_mice;
measurement_name = 'Population vector correlation';
ylims = [0 0.7];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');

clc;
disp('Testing for the effects of time and experience on')
disp('representational drift in population vector correlations')
disp('using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure 3E - the effects of intermediate experience in environment B on ensemble rate correlation values in environment A
figure_name = 'Figure 3E';
stability_measurement = ensemble_rate_context_specificity_all_mice;
measurement_name = 'Ensemble rate correlation';
ylims = [0 0.7];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'context');

clc;
disp('Testing for the effects of intermediate experience in')
disp('environment B on ensemble rate correlation values in')
disp('environment A using two-tailed Wilcoxon signed rank test:')
disp(' ')
disp(statistics_table)

%% Figure 3F - the effects of intermediate experience in environment B on tuning curve correlation values in environment A
figure_name = 'Figure 3F';
stability_measurement = tuning_curve_context_specificity_all_mice;
measurement_name = 'Tuning curve correlation';
ylims = [0 0.7];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'context');

clc;
disp('Testing for the effects of intermediate experience in')
disp('environment B on tuning curve correlation values in')
disp('environment A using two-tailed Wilcoxon signed rank test:')
disp(' ')
disp(statistics_table)

%% Figure 3G - the effects of intermediate experience in environment B on population vector correlation values in environment A
figure_name = 'Figure 3G';
stability_measurement = pv_context_specificity_all_mice;
measurement_name = 'Population vector correlation';
ylims = [0 0.7];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'context');

clc;
disp('Testing for the effects of intermediate experience in')
disp('environment B on population vector correlation values in')
disp('environment A using two-tailed Wilcoxon signed rank test:')
disp(' ')
disp(statistics_table)

%% Figure 4A - Tuning curve stability across days of few example cells in each environment

example_cells_table = array2table([1,1,147; 1,6,169; 1,1,1255; 1,2,1235; 1,1,319; 1,6,256;...
    2 2 901; 2,7,271; 2,5,1058; 2,2,286; 2,6,1136; 2,4,295],...
    'VariableNames',{'environment_id','mouse_id','cell_id'});


figure('units','normalized','position',[0.25 0.3 0.45 0.45])
for cell_ind = 1:size(example_cells_table,1)
    
    env_id = example_cells_table{cell_ind,1};
    mouse_id = example_cells_table{cell_ind,2};
    cell_id = example_cells_table{cell_ind,3};
    
    current_cell = squeeze(population_vectors_all_mice{mouse_id,env_id}(cell_id,:,:))';
    current_cell_norm_tuning = current_cell./max(current_cell,[],2);
    
    number_of_days = size(current_cell,1);
    
    subplot(2,6,cell_ind)
    hold on
    for day_ind = 1:number_of_days
        
        x = [1:20]';
        y = [current_cell_norm_tuning(day_ind,:)+number_of_days-1*(day_ind-1)]';
        fill([x;flipud(x)],[y;ones(size(y))*min(y)],[0.9 0.9 0.9],'linestyle','none');
        plot(y,'color',[0.4 0.4 0.4],'linewidth',1.5)
    end
    ylim([0.975 number_of_days+1])
    xlim([0.85 20.15])
    set(gca,'xtick',[])
    if cell_ind == 1 || cell_ind == 7
        set(gca,'ytick',1.5:number_of_days+1,'yticklabels',fliplr(day_ind_per_env{env_id}))
        ylabel('Day')
    else
        set(gca,'ytick',[])
    end
    if cell_ind == 9
        xlabel('Position on the track')
    end
    
    title(['Cell ',num2str(cell_ind),':'],'FontSize',12,'Fontweight','normal');
    
end
suptitle('Figure 4A')

%% Figure 4B - Primary peak displacement as a function of elpased time for each environment

figure('units','normalized','position',[0.15 0.25 0.5 0.4])
axes('pos',[0 1 1 1],'visible','off');
text(0.425,-0.05,'Figure 4B','fontsize',24);
jet_colormap = colormap(jet);
for env_ind =  1:2
    
    current_env_peak_shift = frac_cells_per_pos_shift_elapsed_time_all_mice(:,env_ind);
    current_env_peak_shift_shuffled = frac_cells_per_pos_shift_elapsed_time_all_mice_shuffled(:,env_ind);
    
    current_mouse_peak_shift = [];
    current_mouse_peak_shift_shuffled = [];
    for mouse_ind = 1:number_of_mice
        current_mouse_peak_shift(:,:,mouse_ind) = current_env_peak_shift{mouse_ind};
        current_mouse_peak_shift_shuffled(:,:,mouse_ind) = current_env_peak_shift_shuffled{mouse_ind};
    end
    
    mean_peak_shift_across_mice = mean(current_mouse_peak_shift,3,'omitnan');
    mean_peak_shift_across_mice_shuffled = mean(current_mouse_peak_shift_shuffled,3,'omitnan');
    
    
    color_factor = 1;
    if env_ind == 1
        color_factor = 0.5;
    end
    
    plt = [];
    lgd_labels = {};
    
    sub1 = axes();
    set(sub1,'units','normalized','position',[0.175 0.15 0.3 0.65]+(env_ind-1)*[0.4 0 0 0])
    hold on
    for interval_ind = 1:size(mean_peak_shift_across_mice,1)
        elapsed_time_in_days = elapsed_days_per_env{env_ind}(interval_ind);
        current_interval_data = mean_peak_shift_across_mice(interval_ind,:);
        current_interval_shuffled = mean_peak_shift_across_mice_shuffled(interval_ind,:);
        
        data_colors = jet_colormap(1+14*color_factor*(interval_ind-1),:)*0.9;
        shuffle_colors = [0.2 0.2 0.2]+0.15*color_factor*(interval_ind-1);
        
        plt(interval_ind) = plot(current_interval_data,'color',data_colors,'linewidth',1.5);
        plot(current_interval_shuffled,'color',shuffle_colors,'linewidth',1.5);
        
        lgd_labels{interval_ind} = num2str(elapsed_time_in_days);
    end
    
    ylim([0 0.225])
    lgd = legend(plt,lgd_labels);
    legend('boxoff')
    title(lgd,'\Delta Days','FontSize',12,'Fontweight','normal')
    
    
    xlabel('Positional shift (cm)')
    set(gca,'xtick',2:6:39,'xticklabels',-72:24:72)
    title([env_full_names{env_ind},':'],'units','normalized','Position',[0.25 1.01],'FontSize',16,'Fontweight','normal');
    
    if env_ind == 1
        ylabel('Fraction of cells')
    elseif env_ind == 2
        set(gca,'yticklabels',[])
    end
    
end


%% Figure 4C - Fraction of stable cells as a function of elapsed time or number of sessions between visits

figure_name = 'Figure 4C';
stability_measurement = fraction_stable_cells_elapsed_time_all_mice;
measurement_name = 'Fraction of stable cells';
ylims = [0 0.6];

[statistics_table, plot_legends] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');

% visualizing shuffle data
envA_shuffle_vals = cell2mat(fraction_stable_cells_elapsed_time_all_mice_shuffled(:,1));
envB_shuffle_vals = cell2mat(fraction_stable_cells_elapsed_time_all_mice_shuffled(:,2));
mean_shuffle_val = mean([envA_shuffle_vals(:);envB_shuffle_vals(:)],'omitnan');
for subplt = 1:2
    subplot(1,2,subplt)
    hold on
    plot(xlim,[mean_shuffle_val,mean_shuffle_val],'--','color',[0.2 0.2 0.2])
    text(0.835,0.2,'Shuffle','units','normalized','fontsize',10);
end
legend(plot_legends)

clc;
disp('Testing for the effects of time and experience on')
disp('the fraction of cells with stable place fields')
disp('using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure 4D - Absolute positional shift in primary place field as a function of elapsed time or number of sessions between visits

figure_name = 'Figure 4D';
stability_measurement = positional_shift_elapsed_time_all_mice;
measurement_name = 'Absolute positional shift (cm)';
ylims = [0 27];

[statistics_table, plot_legends] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');

% visualizing shuffle data
envA_shuffle_vals = cell2mat(positional_shift_elapsed_time_all_mice_shuffled(:,1));
envB_shuffle_vals = cell2mat(positional_shift_elapsed_time_all_mice_shuffled(:,2));
mean_shuffle_val = mean([envA_shuffle_vals(:);envB_shuffle_vals(:)],'omitnan');
for subplt = 1:2
    subplot(1,2,subplt)
    hold on
    plot(xlim,[mean_shuffle_val,mean_shuffle_val],'--','color',[0.2 0.2 0.2])
    text(0.825,0.85,'Shuffle','units','normalized','fontsize',10);
end
legend(plot_legends)

clc;
disp('Testing for the effects of time and experience on')
disp('drift in absolute positional shift in primary place field')
disp('using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure 4E - Absolute change in number of place fields per place cell as a function of elapsed time or number of sessions between visits

figure_name = 'Figure 4E';
stability_measurement = mean_field_number_diff_elapsed_time_all_mice;
measurement_name = 'Absolute change in number of fields';
ylims = [0 1.15];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');

clc;
disp('Testing for the effects of time and experience on')
disp('drift in absolute change in number of place fields per ')
disp('place cell using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure 4F - Absolute change in primary place field size per place cell as a function of elapsed time or number of sessions between visits

figure_name = 'Figure 4F';
stability_measurement = mean_field_width_diff_elapsed_time_all_mice;
measurement_name = 'Absolute change in field size (cm)';
ylims = [0 8];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');

clc;
disp('Testing for the effects of time and experience on')
disp('drift in absolute change in primary place field size per')
disp('place cell using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)


%% Figure S1A - Mean activity rate across the place cells' population over days

figure_name = 'Figure S1A';
stability_measurement = mean_activity_rate_per_sess_all_mice;
measurement_name = 'Mean activity rate (Hz)';
ylims = [0 0.05];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'controls');

clc;
disp('Testing the stability over days of the global mean activity rate across the')
disp('entire place cells population using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure S1B - Number of place cells recorded in each day

figure_name = 'Figure S1B';
stability_measurement = num_place_cells_each_sess_all_mice;
measurement_name = 'Number of place cells';
ylims = [0 300];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'controls');

clc;
disp('Testing the stability over days of number of place cells recorded on')
disp('each day of the experiment using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure S1C - Fraction of place cells recorded in each day

figure_name = 'Figure S1C';
stability_measurement = fraction_pc_each_sess_all_mice;
measurement_name = 'Fraction of place cells';
ylims = [0 0.7];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'controls');

clc;
disp('Testing the stability over days of fraction of place cells recorded on')
disp('each day of the experiment using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure S1D - Fraction of place cells recorded in each day

figure_name = 'Figure S1D';
stability_measurement = mean_info_each_sess_all_mice;
measurement_name = 'Spatial information (bits/event)';
ylims = [0 5];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'controls');

clc;
disp('Testing the stability over days of the mean spatial information across the')
disp('entire place cells population using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure S1E - Average size of primary place field across the place cells' population over days

figure_name = 'Figure S1E';
stability_measurement = mean_place_field_width_all_mice;
measurement_name = 'Primary place field size (cm)';
ylims = [0 17];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'controls');

clc;
disp('Testing the stability over days of the mean place field size across the')
disp('entire place cells population using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure S1F - Average number of place fields per place cell across the place cells' population over days

figure_name = 'Figure S1F';
stability_measurement = mean_place_field_number_all_mice;
measurement_name = 'Number of fields per place cell';
ylims = [0 2.5];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'controls');

clc;
disp('Testing the stability over days of the mean number of fields per cell across the')
disp('entire place cells population using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure S1G - Distribution of number of place fields per place cell in each environment

% Distribution of place fields per place cell for environment A
fraction_cells_per_field_envA = cell2mat(fraction_cells_per_field_num_all_mice(:,1));
mean_fraction_cells_per_field_envA = mean(fraction_cells_per_field_envA,'omitnan');
ste_fraction_cells_per_field_envA = std(fraction_cells_per_field_envA,'omitnan')./sqrt(number_of_mice);

% Distribution of place fields per place cell for environment B
fraction_cells_per_field_envB = cell2mat(fraction_cells_per_field_num_all_mice(:,2));
mean_fraction_cells_per_field_envB = mean(fraction_cells_per_field_envB,'omitnan');
ste_fraction_cells_per_field_envB = std(fraction_cells_per_field_envB,'omitnan')./sqrt(number_of_mice);

plt = [];
figure('units','normalized','position',[0.3 0.25 0.2 0.4])
hold on
plt(1) = bar(1.15:1:4.25,mean_fraction_cells_per_field_envA,'facecolor',env_colors(1,:),'edgecolor','none','barwidth',0.3);
plt(2) = bar(0.85:1:4.5,mean_fraction_cells_per_field_envB,'facecolor',env_colors(2,:),'edgecolor','none','barwidth',0.3);
errorbar(0.85:1:4.5,mean_fraction_cells_per_field_envB,ste_fraction_cells_per_field_envB,'linestyle','none','color','k','capsize',3)
errorbar(1.15:1:4.25,mean_fraction_cells_per_field_envA,ste_fraction_cells_per_field_envA,'linestyle','none','color','k','capsize',3)
ylabel('Fraction of place cells')
xlabel('Number of place fields')
xlim([0.25 4.75])
ylim([0 1])
set(gca,'xtick',1:5)
legend(plt,{'Environment A','Environment B'})
legend('boxoff')
axis square
title('Figure S1G','FontSize',24,'Fontweight','normal')

% sorting data for statistical testing
fraction_cells_per_field_envA_sorted = fraction_cells_per_field_envA';
fraction_cells_per_field_envB_sorted = fraction_cells_per_field_envB';
fraction_cells_per_field_both_env = [fraction_cells_per_field_envA_sorted(:),fraction_cells_per_field_envB_sorted(:)];

[~,stats] = mackskill(fraction_cells_per_field_both_env,4);

field_num_stats =  [stats.T,stats.df,stats.p];


col_names = {'Chi','df','pvalue'};
statistics_table = array2table(field_num_stats,'VariableNames',col_names);

clc;
disp('Testing the similarity between the distributions of number of place fields ')
disp('per place cell across the two environments using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)


%% Figure S1H - Population vector correlation between environments as a function of time

% similarity between visits of different environments
mean_pv_corr_across_env_all_mice_sorted = mean_pv_corr_across_env_all_mice; % to be removed
mean_pv_corr_across_env_and_mice = mean(mean_pv_corr_across_env_all_mice_sorted,'omitnan');
ste_pv_corr_across_env_and_mice = std(mean_pv_corr_across_env_all_mice_sorted,'omitnan')./sqrt(number_of_mice);

% similiarty between visits of the same environment
pv_corr_elapsed_time_all_mice_envA = cell2mat(pv_corr_elapsed_time_all_mice(:,1));
pv_corr_between_consecutive_sess = pv_corr_elapsed_time_all_mice_envA(:,1);
same_env_stability = mean(pv_corr_between_consecutive_sess,'omitnan');

figure('units','normalized','position',[0.3 0.25 0.2 0.4])
xlim([0 11])
ylim([-0.075 0.7])
hold on
% plot data from different environments
plot(mean_pv_corr_across_env_all_mice_sorted','color',[0.8 0.8 0.8])
errorbar(mean_pv_corr_across_env_and_mice,ste_pv_corr_across_env_and_mice,'-o','color','k','linewidth',1.5,'MarkerSize',4,...
    'MarkerEdgeColor','k','MarkerFaceColor','k');

% plot data from the same environment
plot(xlim,[same_env_stability, same_env_stability],'--','color',[0.2 0.2 0.2])
text(0.6,0.825,'Same environment','units','normalized','fontsize',10);

set(gca,'xtick',1:10,'xticklabels',day_ind_per_env{1})
xlabel('Day')
ylabel({'Population vector correlation';'between different environments'})
axis square
title('Figure S1H','FontSize',24,'Fontweight','normal')

% performing statistical testing
[~,drift_stats] = mackskill(mean_pv_corr_across_env_all_mice_sorted,1);
col_names = {'Chi','df','pvalue'};
statistics_table = array2table([drift_stats.T,drift_stats.df,drift_stats.p],'VariableNames',col_names);

clc;
disp('Testing the stability of the similarity between the representations of ')
disp('different environments over time using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure S1I - Effects of time and experiences on representational drift when including reward sites

% Repeating the ensemble rate correlation analysis (see Figure 3B) while including reward
figure_name = 'Figure S1I (top)';
stability_measurement = rate_corr_elapsed_time_with_edges_all_mice;
measurement_name = 'Ensemble rate correlation';
ylims = [0 0.7];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');

clc;
disp('Testing the effects of time and experience on')
disp('representational drift when including reward sites')
disp('using the nonparametric Skillings-Mack test:')
disp(' ')
disp('Ensemble rate correlation:')
disp(statistics_table)

% Repeating the tuning curve correlation analysis (see Figure 3C) while including reward
figure_name = 'Figure S1I (bottom)';
stability_measurement = tuning_corr_elapsed_time_with_edges_all_mice;
measurement_name = 'Tuning curve correlation';
ylims = [0 0.7];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');

disp(' ')
disp('Tuning curve correlation:')
disp(' ')
disp(statistics_table)


%% Figure S1J - Effects of time and experiences on representational drift when including reward sites

% Repeating the ensemble rate correlation analysis (see Figure 3B)
% while including only the reward sites
figure_name = 'Figure S1J (top)';
stability_measurement = rate_corr_elapsed_time_only_edges_all_mice;
measurement_name = 'Ensemble rate correlation';
ylims = [0 0.7];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');

clc;
disp('Testing the effects of time and experience on')
disp('representational drift while including only reward sites')
disp('using the nonparametric Skillings-Mack test:')
disp(' ')
disp('Ensemble rate correlation:')
disp(statistics_table)

% Repeating the tuning curve correlation analysis (see Figure 3C)
% while including only the reward sites
figure_name = 'Figure S1J (bottom)';
stability_measurement = tuning_corr_elapsed_time_only_edges_all_mice;
measurement_name = 'Tuning curve correlation';
ylims = [0 0.7];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');

disp(' ')
disp('Tuning curve correlation:')
disp(' ')
disp(statistics_table)

%% Figure S1K - Effects of time and experiences on representational drift when including reward sites

% Repeating analysss presented in Figures 3B-3B while controling
% for the number of registered days across environments
figure_name = 'Figure S1K';
stability_measurement = {[ensemble_rate_corr_elapsed_time_all_mice(:,1),rate_corr_elapsed_time_registration_control_all_mice],...
    [tuning_curve_corr_elapsed_time_all_mice(:,1),tuning_corr_elapsed_time_registration_control_all_mice]};
measurement_name = {'Ensemble rate correlation','Tuning curve correlation'};
ylims = [0 0.7];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'registration');

clc;
disp('Testing the effects of time and experience on')
disp('representational drift in ensemble rate and tuning curve ')
disp('correlations after controling for the number of registered days ')
disp('across environment using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure S1L -  Relationship between absolute activity rate difference and tuning curve stability at the single cell resolution

figure_name = 'Figure S1L';
stability_measurement = tuning_rate_measurements_all_mice;
measurement_name = 'Absolute activity rate difference (Hz)';
input_data = {1,[0 0.7]};

time_vs_exp_main_visual(stability_measurement,input_data,measurement_name,figure_name,'independence');

%% Figure S1M -  Relationship between absolute activity rate difference score and tuning curve stability at the single cell resolution

figure_name = 'Figure S1M';
stability_measurement = tuning_rate_measurements_all_mice;
measurement_name = 'Absolute activity rate difference score';
input_data = {2,[0 1.1]};

time_vs_exp_main_visual(stability_measurement,input_data,measurement_name,figure_name,'independence');

%% Figure S1N -  Relationship between mean activity rate and tuning curve stability at the single cell resolution

figure_name = 'Figure S1N';
stability_measurement = tuning_rate_measurements_all_mice;
measurement_name = 'Mean activity rate (Hz)';
input_data = {3,[0 0.8]};

time_vs_exp_main_visual(stability_measurement,input_data,measurement_name,figure_name,'independence');


%% Figure S1O -  Relationship between mean activity rate and tuning curve stability at the single cell resolution

% averaging the correlations of each mouse across the two environments
avg_tuning_rate_relationship = [];
for mouse_ind = 1:number_of_mice
    current_mouse_both_env = cell2mat(tuning_rate_relationship_all_mice(mouse_ind,:)');
    avg_current_mouse = mean(current_mouse_both_env,'omitnan');
    
    avg_tuning_rate_relationship(mouse_ind,:) = avg_current_mouse;
end


figure('units','normalized','position',[0.3 0.25 0.2 0.35])
figure_boxplot(avg_tuning_rate_relationship);
ylim([0 1])
set(gca,'xtick',1:3,'xticklabels',{'Abs. rate difference','Abs. rate diff. score','Mean activity rate'})
xtickangle(30)
xlim([0 4])
ylabel({'Explained variance by';'tuning curve correlation (R^2)'})
axis square
title('Figure S1O','FontSize',24,'Fontweight','normal')

%% Figure S2A - Population vector correlation across sessions for each environment at the single animal resolution

figure('units','normalized','position',[0.2075 0 0.1175 1])
axes('pos',[0 1 1 1],'visible','off');
text(0.175,-0.035,'Figure S2A','fontsize',24);

mouse_group_id = [1,1,1,2,2,2,1,1];
environment_types = {'Env. A (Straight)','Env. B (L-shaped)'; 'Env. A (L-shaped)','Env. B (Straight)'};

for env_ind = 1:2
    current_env_days = day_ind_per_env{env_ind};
    current_env_pv_corr_mat = cell2mat(pv_corr_mat_all_mice(:,env_ind,:));
    for mouse_ind = 1:number_of_mice
        current_mouse_pv_corr_mat = current_env_pv_corr_mat(:,:,mouse_ind);
        
        subplot(number_of_mice,2,env_ind+2*(mouse_ind-1))
        imagesc(current_mouse_pv_corr_mat,'AlphaData',~isnan(current_mouse_pv_corr_mat),[-0.1 1])
        hold on
        for line = 1:(size(current_mouse_pv_corr_mat,1)./20)-1
            plot([20.5 20.5]+20*(line-1),ylim,'color','k')
            plot(xlim,[20.5 20.5]+20*(line-1),'color','k')
        end
        set(gca,'xtick',[],'ytick',[],'color',[1 1 1])
        title(environment_types{mouse_group_id(mouse_ind),env_ind},'FontSize',8,'Fontweight','normal')
        if env_ind == 1
            ylabel('Day','FontSize',8)
        end
        
        if mouse_ind == 8
            xlabel('Day','FontSize',8)
        end
    
        
    end
end
colormap(jet)

%% Figure S2B - Effects of time and experience on ensemble rate correlations at the single animal resolution
figure_name = 'Figure S2B';
stability_measurement = ensemble_rate_corr_elapsed_time_all_mice;
measurement_name = 'Ensemble rate corr.';
input_data = [0 0.7];

figure('units','normalized','position',[0.35 0 0.125 1])
time_vs_exp_main_visual(stability_measurement,input_data,measurement_name,figure_name,'individual');

%% Figure S2C - Effects of time and experience on tuning curve correlations at the single animal resolution
figure_name = 'Figure S2C';
stability_measurement = tuning_curve_corr_elapsed_time_all_mice;
measurement_name = 'Tuning curve corr.';
input_data = [0 0.7];

figure('units','normalized','position',[0.5 0 0.125 1])
time_vs_exp_main_visual(stability_measurement,input_data,measurement_name,figure_name,'individual');

%% Figure S2D - Effects of time and experience on tuning curve correlations at the single animal resolution
figure_name = 'Figure S2D';
stability_measurement = pv_corr_elapsed_time_all_mice;
measurement_name = 'PV correlation';
input_data = [0 0.7];

figure('units','normalized','position',[0.65 0 0.125 1])
time_vs_exp_main_visual(stability_measurement,input_data,measurement_name,figure_name,'individual');

%% Figure S3 - Quantifying the relative contribution of changes in activity rates and spatial tuning to representational drift using simulation

% IMPORTANT NOTE:
% It is required to run this section in order to plot Figures S3A-S3F.
% To replicate the results in the published manuscript it is required to
% run exactly 1000 realizations of the simulation with the provided rng seed
% using "drift_simulation" function. This process takes approximately 4:45 hours. 

% define simulation parameters
simulation_parameters = [];
simulation_parameters.number_of_cells = 1000;
simulation_parameters.tuning_smoothing_sigma = 2;
simulation_parameters.tuning_change_sigma = 2;
simulation_parameters.activity_rate_mu = log(0.025);
simulation_parameters.activity_rate_sigma = 0.75;
simulation_parameters.activity_sigma_factor = 3;
simulation_parameters.prob_tuning_change = 0.6;
simulation_parameters.prob_rate_change = 0.75;
simulation_parameters.fraction_of_place_cells = 0.6;
simulation_parameters.number_of_bins = 20;
simulation_parameters.number_of_days = 10;
simulation_parameters.number_of_simulations = 1000;

[full_drift_model,tuning_drift_model,rate_drift_model] = drift_simulation(simulation_parameters,state);


%% Figure S3A - Neuronal responses of six example simulated place cells
cells_across_all_simu  = cell2mat(full_drift_model.population_vectors_each_day');
example_cells = [1813,4749,7112,7946,23277,25143];
number_of_days = simulation_parameters.number_of_days;

figure('units','normalized','position',[0.3 0.3 0.2 0.4])
for cell_ind = 1:length(example_cells)
    
    current_cell =  squeeze(cells_across_all_simu(example_cells(cell_ind),:,:))';
    current_cell_norm = current_cell./max(current_cell,[],2);
    current_cell_norm(isnan(current_cell_norm)) = 0;
    
    
    subplot(2,3,cell_ind)
    for day = 1:10
        hold on
        x = [1:20]';
        y = [current_cell_norm(day,:)+10-1*(day-1)]'-min(current_cell_norm(day,:));
        fill([x;flipud(x)],[y;[ones(size(y))*min(y)]-0.1],[0.9 0.9 0.9],'linestyle','none');
        
        
        plot(y,'color',[0.4 0.4 0.4],'linewidth',1.5)
    end
    

    set(gca,'xtick',[])
    if cell_ind == 1 || cell_ind == 4
        set(gca,'ytick',1.5:number_of_days+1,'yticklabels',fliplr(1:number_of_days))
        ylabel('Day')
    else
        set(gca,'ytick',[])
    end
    if cell_ind == 5
        xlabel('Position on the track')
    end
        ylim([1 number_of_days+1])
    xlim([0.85 20.15])
    title(['Cell ',num2str(cell_ind),':'],'FontSize',12,'Fontweight','normal');
    
end
suptitle('Figure S3A')


%% Figure S3B - Spatial tuning curves of the same simulated place cells across six days

% IMPORTANT NOTE%
% Please note that it is required to save the plots as PDF of highest resolution
% in order to obtain the same visuals as in the published manuscript.

example_siulation_ind = 1000;
number_of_cells = simulation_parameters.number_of_cells;
example_simulation_population_vectors = full_drift_model.population_vectors_each_day{example_siulation_ind};
example_simulation_place_cells_ind = full_drift_model.PC_ind_each_day(example_siulation_ind,:);
reference_day_list = 1:2:10;
day_ind_list = 1:2:10;
sub_plt_ind = 1;

figure('units','normalized','position',[0.3 0.15 0.2 0.55])
axes('pos',[0 1 1 1],'visible','off');
text(0.375,-0.035,'Figure S3B','fontsize',16);
for ref_day_ind = 1:length(reference_day_list)
    reference_day = reference_day_list(ref_day_ind);
    
    reference_day_pc_ind = example_simulation_place_cells_ind{reference_day};
    reference_day_population_vectors = example_simulation_population_vectors(reference_day_pc_ind,:,reference_day);
    [~,ref_day_peak_pos] = max(reference_day_population_vectors,[],2);
    [~,ref_day_sorted_cells_ind] = sort(ref_day_peak_pos);
   
    for day_ind = 1:length(day_ind_list)
        current_day_ind = day_ind_list(day_ind);
        
        current_day_pc_ind = example_simulation_place_cells_ind{current_day_ind};
        current_day_non_pc_ind = ismember(1:number_of_cells,current_day_pc_ind)';
        current_day_tuning_curves = example_simulation_population_vectors(:,:,current_day_ind).*current_day_non_pc_ind;
        current_day_tuning_curves_pc = current_day_tuning_curves(reference_day_pc_ind,:);
        
        current_day_peak_activity = max(current_day_tuning_curves_pc,[],2);
        current_day_norm_tuning_curves = current_day_tuning_curves_pc./current_day_peak_activity;
        current_day_sorted_tuning_curves = current_day_norm_tuning_curves(ref_day_sorted_cells_ind,:);
        
        subplot(length(reference_day_list),length(day_ind_list),sub_plt_ind)
        imagesc(current_day_sorted_tuning_curves)
        
 
        if ref_day_ind == 5
            xlabel(['Day ',num2str(current_day_ind)])
        end
        if day_ind == 1
            ylabel(['Day ',num2str(reference_day)])
        end
        
        if ref_day_ind == 1 && day_ind == 5
            title({'Position'},'units','normalized','Position',[0.5 1.01],'FontSize',10,'Fontweight','normal');
            set(gca,'YAxisLocation','right')
            ylabel('Cell order','FontSize',10)
        end
        set(gca,'xtick', [],'ytick', []);
        
        sub_plt_ind = sub_plt_ind + 1;
    end
end
jet_colormap = colormap(jet);
colormap(jet_colormap*0.9)

cb = colorbar();
cb.Position = [0.925 0.107 0.03 0.4725];
cb.Ticks = [0 1];
ylabel(cb,'Normalized activity rate','FontSize',10,'Rotation',270);



%% Figure S3C - Population vector correlation between the representations of different locations across 10 days.

pv_corr_mat_all_simulations = full_drift_model.pv_corr_mat;
avg_pv_corr_mat_across_simulations = mean(pv_corr_mat_all_simulations,3,'omitnan');
number_of_bins = simulation_parameters.number_of_bins;

figure('units','normalized','position',[0.25 0.2 0.25 0.35])
axes('pos',[0 1 1 1],'visible','off');
text(0.375,-0.05,'Figure S3C','fontsize',24);
text(0.69,-0.15,'position','fontsize',9);

sub1 = axes();
set(sub1,'units','normalized','position',[0.25 0.125 0.525 0.7])
imagesc(avg_pv_corr_mat_across_simulations,[0.1 1])
hold on
for line = 1:(size(avg_pv_corr_mat_across_simulations,1)./number_of_bins)-1
    plot([20.5 20.5]+number_of_bins*(line-1),ylim,'color','k')
    plot(xlim,[20.5 20.5]+number_of_bins*(line-1),'color','k')
end
set(gca,'xtick',10.5:20:size(avg_pv_corr_mat_across_simulations,1),'xticklabel',1:number_of_days,...
    'ytick',10.5:20:size(avg_pv_corr_mat_across_simulations,1),'yticklabel',1:number_of_days)
xlabel('Day','FontSize',12)
ylabel('Day','FontSize',12)
title(['Full drift model:'],'units','normalized','Position',[0.22 1.01],'FontSize',12,'Fontweight','normal');
cb = colorbar();
cb.Position = [0.8 0.125 0.0225 0.35];
cb.Ticks = [0.1 1];
ylabel(cb,'PV correlation','FontSize',12,'Rotation',270);
colormap(jet)

%% Figure S3D - Stability measures as a function of elapsed time for the three different model

measurment_list = {full_drift_model.ensemble_rate_corr_elapsed_time,full_drift_model.tuning_curve_corr_elapsed_time,full_drift_model.pv_corr_elapsed_time;...
    tuning_drift_model.ensemble_rate_corr_elapsed_time,tuning_drift_model.tuning_curve_corr_elapsed_time,tuning_drift_model.pv_corr_elapsed_time;...
    rate_drift_model.ensemble_rate_corr_elapsed_time,rate_drift_model.tuning_curve_corr_elapsed_time,rate_drift_model.pv_corr_elapsed_time};

measurment_names = {'Ensemble rate corr.','Tuning curve corr.','PV  correlation'};
model_names = {'Full drift model','Tuning drift model','Rate drift model'};

figure('units','normalized','position',[0.2 0.1 0.375 0.65])

for model_ind = 1:3
    for measurment_ind = 1:3
        
        
        current_measurment = measurment_list{model_ind,measurment_ind};
        mean_stability = nanmean(current_measurment);
        ste_stability = nanstd(current_measurment);
        
        
        subplot(3,3,model_ind+3*(measurment_ind-1))
        hold on
        
        errorbar(mean_stability,ste_stability,'-o','color','k','linewidth',1.5,'MarkerSize',4,...
            'MarkerEdgeColor','k','MarkerFaceColor','k');
        ylim([0 0.7])
        xlim([0 10])
        axis square
        
        if measurment_ind == 1
            title([model_names{measurment_ind},':'],'units','normalized','FontSize',12,'Fontweight','normal');
        end
        
        if measurment_ind == 3
            xlabel('Days difference','FontSize',10)
            set(gca,'xtick',1:2:9)
        else
            set(gca,'xticklabels',[])
        end
        
        if model_ind == 1
            ylabel(measurment_names{measurment_ind},'FontSize',10)
            set(gca,'ytick',0:0.1:0.7,'FontSize',10)
        else
            set(gca,'yticklabels',[])
        end
        
        
    end
end

suptitle('Figure S3D')

%% Figure S3E - Population vector correlation as function of elapsed time for the three different model

% load required colormap
load(['Y:\personal\nitzang\time_v_exp\neuron_revisions_2022\scripts\magma_colormap.mat'])


% Population vector correlation for the three models
measurment_list = {full_drift_model.pv_corr_elapsed_time,...
    tuning_drift_model.pv_corr_elapsed_time,...
    rate_drift_model.pv_corr_elapsed_time};


plt = [];
figure('units','normalized','position',[0.3 0.25 0.2 0.4])
for model_ind = 1:3
    current_measurment = measurment_list{model_ind};
    mean_stability = nanmean(current_measurment);
    ste_stability = nanstd(current_measurment);
    
    hold on
    plt(model_ind) = errorbar(mean_stability,ste_stability,'-o','color',magma_colormap(1+110*(model_ind-1),:),'linewidth',2,'MarkerSize',5,...
        'MarkerEdgeColor',magma_colormap(1+110*(model_ind-1),:),'MarkerFaceColor',magma_colormap(1+110*(model_ind-1),:));
    xlim([0 10])
    ylim([0 0.7])
    set(gca,'xtick',1:2:9)
    axis square
    xlabel('Days difference')
    
    ylabel('Population vector correlation')
end
title('Figure S3E','FontSize',24,'Fontweight','normal')
legend(plt,{'Full drift','Tuning drift','Rate drift'},'Location','southwest')
legend('boxoff')
axis square

%% Figure S3F - Population vector correlation as function of elapsed time for the three different model

% Decoding performance for the three models
model_list = {full_drift_model.decoder_accuracy_elapsed_time,...
    tuning_drift_model.decoder_accuracy_elapsed_time,...
    rate_drift_model.decoder_accuracy_elapsed_time};

plt = [];
figure('units','normalized','position',[0.3 0.25 0.2 0.4])
ylim([0 100])
xlim([0 10])
hold on
plot(xlim,[100 100]./number_of_bins,'--','color',[0.2 0.2 0.2],'linewidth',1.5)
text(0.35,9,'Chance','fontsize',10);
for model_ind = 1:3
    current_model = model_list{model_ind}*100;
    mean_acc = nanmean(current_model);
    ste_acc = nanstd(current_model);
    
    
    plt(model_ind) = errorbar(mean_acc,ste_acc,'-o','color',magma_colormap(1+110*(model_ind-1),:),'linewidth',2,'MarkerSize',5,...
        'MarkerEdgeColor',magma_colormap(1+110*(model_ind-1),:),'MarkerFaceColor',magma_colormap(1+110*(model_ind-1),:));
    set(gca,'xtick',1:2:9,'ytick',0:20:100)
    axis square
    xlabel('Days difference')
    ylabel('% of correct classifications')
end
title('Figure S3F','FontSize',24,'Fontweight','normal')

legend(plt,{'Full drift','Tuning drift','Rate drift'},'Location','northwest')
legend('boxoff')


%% Figure S4A - Position along the linear track in two example sessions in Environment B for a single example mouse

figure_name = 'Figure S4A';
env_ind = 2; % Environment B
example_mouse_ind = 5;
velocity_threshold = 2; % cm/ses
bins_to_use = [2:22]; % exluding reward zones bins 1-2 and bins 23-24
example_sess_list = [3,6]; % Day 9 and day 21
frame_rate = 10; % Hz
xticks_in_frames = [1200:1200:12000]; 
xtick_in_seconds = xticks_in_frames./frame_rate;
xtick_in_minutes = xtick_in_seconds./60;

example_mouse_location_per_sess = location_per_sess_all_mice{example_mouse_ind,env_ind};
example_mouse_velocity_per_sess = velocity_per_sess_all_mice{example_mouse_ind,env_ind};

figure('units','normalized','position',[0.1 0.25 0.7 0.5])
axes('pos',[0 1 1 1],'visible','off');
text(0.45,-0.035,figure_name,'fontsize',20);

for sess_ind = 1:length(example_sess_list)
    current_sess_ind = example_sess_list(sess_ind);
    current_day_ind = day_ind_per_env{env_ind}(current_sess_ind);
    
    current_sess_locations = example_mouse_location_per_sess{current_sess_ind};
    current_sess_velocity = example_mouse_velocity_per_sess{current_sess_ind};
  
    if sess_ind == 1
       xlims = [-250 length(current_sess_locations)+250];
    end
    
    subplot(2,1,sess_ind)
    hold on
    plot(current_sess_locations,'color',[0.5 0.5 0.5])
    plot([7000 7000],ylim,'k','linewidth',1.25)
    plot([7625 7625],ylim,'k','linewidth',1.25)
    plot([7000 7625],[max(ylim) max(ylim)],'k','linewidth',1.25)
    
    ylim([0.5 29.5])
    xlim(xlims)
    ylabel('Position (bin)','FontSize',14)
    set(gca,'xtick',xticks_in_frames,'xticklabels',xtick_in_minutes,'ytick',5:5:20)
    box off
    
    if sess_ind == 1
        set(gca,'xtick',[])
    elseif sess_ind == 2
        xlabel('Time (min)','FontSize',14)
    end
     text(0.015, 0.9,['Session ',num2str(current_sess_ind),' (day ',num2str(current_day_ind),'):'],...
         'units','normalized','FontSize',14);
   
end

%% Figure S4B - Magnification of a few representative track traversals from the same example sessions presented in Figure S4A

figure_name = 'Figure S4B';
env_ind = 2; % Environment B
example_mouse_ind = 5;
velocity_threshold = 2; % cm/ses
bins_to_use = [3:22]; % exluding reward zones bins 1-2 and bins 23-24
example_sess_list = [3,6]; % Day 9 and day 21
frame_rate = 10; % Hz
xticks_in_frames = [7050:150:7625]; 
xtick_in_seconds = xticks_in_frames./frame_rate;
xtick_in_minutes = xtick_in_seconds./60;

xlims = [7000 7625];

example_mouse_location_per_sess = location_per_sess_all_mice{example_mouse_ind,env_ind};
example_mouse_velocity_per_sess = velocity_per_sess_all_mice{example_mouse_ind,env_ind};

figure('units','normalized','position',[0.3 0.25 0.25 0.4])
axes('pos',[0 1 1 1],'visible','off');
text(0.625,-0.035,figure_name,'fontsize',20);
text(0.125,-0.02,'Right running direction','fontsize',10,'color',magma_colormap(180,:));
text(0.125,-0.055,'Left running direction','fontsize',10,'color',magma_colormap(115,:));
text(0.95,-0.125,'Session 3 (day 9)','fontsize',12,'rotation',270);
text(0.95,-0.55,'Session 6 (day 21)','fontsize',12,'rotation',270);

for sess_ind = 1:length(example_sess_list)
    current_sess_ind = example_sess_list(sess_ind);
    current_day_ind = day_ind_per_env{env_ind}(current_sess_ind);
    
    current_sess_locations = example_mouse_location_per_sess{current_sess_ind};
    current_sess_velocity = example_mouse_velocity_per_sess{current_sess_ind};
    current_sess_absolute_velocity = abs(current_sess_velocity);
    
    in_track_frames = ismember(current_sess_locations,bins_to_use);
    right_running_direction_ind = find(current_sess_velocity >= velocity_threshold & in_track_frames);
    left_running_direction_ind = find(current_sess_velocity <= -velocity_threshold & in_track_frames);
    
    current_sess_loc_right_run_dir = current_sess_locations(right_running_direction_ind);
    current_sess_vel_right_run_dir = current_sess_velocity(right_running_direction_ind);
    
    current_sess_loc_left_run_dir = current_sess_locations(left_running_direction_ind);
    current_sess_vel_left_run_dir = current_sess_velocity(left_running_direction_ind);
    
    subplot(4,1,1+2*(sess_ind-1))
    hold on
    plot(current_sess_locations,'color',[0.5 0.5 0.5])
    scatter(right_running_direction_ind,current_sess_loc_right_run_dir,10,magma_colormap(180,:),'filled')
    scatter(left_running_direction_ind,current_sess_loc_left_run_dir,10,magma_colormap(115,:),'filled')
    box off
    ylim([0.5 24.5])
    xlim(xlims)
    ylabel({'Position';'(bin)'})
    set(gca,'xticklabels',[])
    
    
    subplot(4,1,2+2*(sess_ind-1))
    hold on
    plot(current_sess_absolute_velocity ,'color',[0.5 0.5 0.5])
    box off
    ylim([0 75])
     xlim(xlims)
    ylabel({'Velocity';'cm/sec)'})
    set(gca,'ytick',[0:25:75])

    if sess_ind == 2
     set(gca,'xtick',xticks_in_frames,'xticklabels',xtick_in_minutes)
     xlabel('Time (min)')
    elseif sess_ind == 1
        set(gca,'xtick',[])
    end
end

%% Figure S4C - Example location occupancy profile and velocity profile along the track across days for each of the two environments

figure_name = 'Figure S4C';
example_mouse_ind = 5;

example_mouse_occupancy_profile = occupancy_profile_per_sess_all_mice(example_mouse_ind,:);
example_mouse_velocity_profile = velocity_profile_per_sess_all_mice(example_mouse_ind,:);

figure('units','normalized','position',[0.2 0.25 0.4 0.3])
axes('pos',[0 1 1 1],'visible','off');
text(0.425,-0.03,figure_name,'fontsize',18);

colormap(magma_colormap)
cb = colorbar();
cb.Position = [0.925 0.125 0.015 0.8];
cb.Ticks = [0 1];
cb.YTickLabel = {'First','Last'};
ylabel(cb,'Day of experiment','FontSize',12,'Rotation',270);

for env_ind = 1:2
	current_env_occupancy_profile = example_mouse_occupancy_profile{env_ind};
    current_env_velocity_profile = example_mouse_velocity_profile{env_ind};

    number_of_sessions = size(current_env_occupancy_profile,1);
    for sess_ind = 1:number_of_sessions

        current_sess_occupancy_profile = current_env_occupancy_profile(sess_ind,:);
        current_sess_velocity_profile = abs(squeeze(current_env_velocity_profile(sess_ind,:,1)));
        
        % fliping the velocity profiles to have the same direction for visulization purposes
        current_sess_velocity_profile = [current_env_velocity_profile(sess_ind,:,1);...
                                         fliplr(current_env_velocity_profile(sess_ind,:,2))];
                                     
        mean_velocity_profile_across_run_dir = mean(current_sess_velocity_profile,'omitnan');
            
        if env_ind == 1 % Environment A (straight track)
            current_sess_color = magma_colormap(216-22*(sess_ind-1),:);
            % averaging the bins 12 and 13 in the L-shaped track
            % (connection point between arms) for visualization purposes
            current_sess_occupancy_profile = [current_sess_occupancy_profile(1:11),...
                                            sum(current_sess_occupancy_profile(12:13),'omitnan'),...
                                            current_sess_occupancy_profile(14:end)];

             mean_velocity_profile_across_run_dir = [mean_velocity_profile_across_run_dir(1:11),...
                                                    mean(mean_velocity_profile_across_run_dir(12:13),'omitnan'),...
                                                    mean_velocity_profile_across_run_dir(14:end)];
            xticks = [1 12 23];
        elseif env_ind == 2 % Environment B (L-shaped track)
            current_sess_color = magma_colormap(216-40*(sess_ind-1),:);
            xticks = [1 12.5 24];
          
            end
        
        subplot(2,2,1+(env_ind-1))
        hold on
        plot(current_sess_occupancy_profile,'color',current_sess_color,'linewidth',2);
        set(gca,'xtick',xticks,'xticklabels',[],'ytick',0:0.05:0.2)
        ylim([0 0.2])
        
        if env_ind == 1
        ylabel('Fraction of time')
        else
            set(gca,'yticklabels',[])
        end
        
        subplot(2,2,3+(env_ind-1))
        hold on
        plot(mean_velocity_profile_across_run_dir,'color',current_sess_color,'linewidth',2);
        ylim([0 50])
        set(gca,'xtick',xticks,'xticklabels',[0 48 96],'ytick',0:10:50)
        xlabel('Position (cm)')
        
        if env_ind == 1
            ylabel('Velocity (cm/sec)')
        else
            set(gca,'yticklabels',[])
        end
    end
end

%% Figure S4D - Occupancy profile similarity as a function of elapsed time or number of sessions between visits
figure_name = 'Figure S4D';
stability_measurement = occupancy_profile_elapsed_time_all_mice;
measurement_name = 'Occupancy profile similarity (r)';
ylims = [0 1];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');

clc;
disp('Testing for the effects of time and experience on')
disp('occupancy profile similarity using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure S4E - Velocity profile similarity as a function of elapsed time or number of sessions between visits
figure_name = 'Figure S4E';
stability_measurement = velocity_profile_elapsed_time_all_mice;
measurement_name = 'Velocity profile similarity (r)';
ylims = [0 1];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');

clc;
disp('Testing for the effects of time and experience on')
disp('velocity profile similarity using the nonparametric Skillings-Mack test:')
disp(' ')
disp(statistics_table)

%% Figure S4E -  Effects of time and experience on representational drift after controling for differences in number of traversals between visits

figure_name = 'Figure S4F (top)';
stability_measurement = rate_corr_elapsed_time_all_mice_behavioral_control;
measurement_name = 'Ensemble rate correlation';
ylims = [0 0.7];

[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');

clc;
disp('Testing for the effects of time and experience on representational drift')
disp('after controling for differences in number of traversals between visits')
disp('using the nonparametric Skillings-Mack test:')
disp(' ')
disp('Ensemble rate correlation:')
disp(statistics_table)

figure_name = 'Figure S4F (middle)';
stability_measurement = tuning_corr_elapsed_time_all_mice_behavioral_control;
measurement_name = 'Tuning curve correlation';
[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');


disp(' ')
disp('Tuning curve correlation:')
disp(statistics_table)

figure_name = 'Figure S4F (bottom)';
stability_measurement = pv_corr_elapsed_time_all_mice_behavioral_control;
measurement_name = 'Population vector correlation';
[statistics_table] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,'main');


disp(' ')
disp('Population vector correlation:')
disp(statistics_table)


%% Figure S4G - Linear relationship between neuronal code stability and behavioral variability

measurement_names = {'ensemble rate correlation','tuning curve correlation','population vector correlation'};
xtick_labels = {'Time','Experience','Running speed','Number of traversals','Occupancy profile','Velocity profile'};


figure('units','normalized','position',[0.1 0.25 0.75 0.375])
axes('pos',[0 1 1 1],'visible','off');
text(0.475,-0.03,'Figure S4G','fontsize',20);

for measurement_ind = 1:length(measurement_names)
    
    stability_behavior_relationship = [];
    for mouse_ind = 1:number_of_mice
        current_mouse_regression = [behavior_regression_all_mice{mouse_ind,1};behavior_regression_all_mice{mouse_ind,2}];
        stability_behavior_relationship(mouse_ind,:) = corr(current_mouse_regression{:,measurement_ind},current_mouse_regression{:,4:end},'rows','pairwise');
    end
    
    stability_behavior_exaplained_variance = stability_behavior_relationship.^2;

    subplot(1,length(measurement_names),measurement_ind)
    figure_boxplot(stability_behavior_exaplained_variance);
    
    ylabel({'Variability explained by';[measurement_names{measurement_ind}, ' (R^2)']})
    ylim([0 1])
    xlim([-0.75 7.75])
    set(gca,'xtick',1:6,'xticklabels',xtick_labels)
    xtickangle(45)
end
    
