function [statistics_table,plt] = time_vs_exp_main_visual(stability_measurement,ylims,measurement_name,figure_name,flag)

day_ind_per_env = {2:2:20,1:4:21};
env_colors = [10 105 170; 190 25 50]./255;
elapsed_days_per_env = {2:2:18;4:4:20};
env_full_names = {'Environment A','Environment B'};
number_of_mice = size(stability_measurement,1);

if strcmp(flag,'main')
    
    figure('units','normalized','position',[0.25 0.25 0.4 0.4])
    axes('pos',[0 1 1 1],'visible','off');
    text(0.425,-0.05,figure_name,'fontsize',24);
    
    plt = [];
    for env_ind = 1:2
        elapsed_days = elapsed_days_per_env{env_ind};
        
        current_env_stability = cell2mat(stability_measurement(:,env_ind));
        number_of_mice = sum(~isnan(current_env_stability));
        
        mean_stability_across_mice = mean(current_env_stability,'omitnan');
        ste_stability_across_mice = std(current_env_stability,'omitnan')./sqrt(number_of_mice);
        
        subplot(1,2,1)
        hold on
        plt(env_ind) = errorbar(elapsed_days,mean_stability_across_mice, ste_stability_across_mice,'-o','color',env_colors(env_ind,:),'linewidth',2,'MarkerSize',4,...
            'MarkerEdgeColor',env_colors(env_ind,:),'MarkerFaceColor',env_colors(env_ind,:));
        ylim(ylims)
        xlim([0 21])
        axis square
        ylabel(measurement_name)
        set(gca,'xtick',2:2:20)
        if env_ind == 2
            legend(plt,{'Environment A','Environment B'},'Location','southwest')
            legend('boxoff')
        end
        xlabel('Days difference')
        
        subplot(1,2,2)
        hold on
        errorbar(mean_stability_across_mice, ste_stability_across_mice,'-o','color',env_colors(env_ind,:),'linewidth',2,'MarkerSize',4,...
            'MarkerEdgeColor',env_colors(env_ind,:),'MarkerFaceColor',env_colors(env_ind,:));
        ylim(ylims)
        xlim([0 10])
        axis square
        xlabel('Sessions difference')
        set(gca,'xtick',1:9,'yticklabels',[])
    end
    
    
    % sort data for statistical testing
    stability_measurement_envA = cell2mat(stability_measurement(:,1));
    stability_measurement_envB = cell2mat(stability_measurement(:,2));
    
    [~,drift_stats_envA] = mackskill(stability_measurement_envA,1);
    [~,drift_stats_envB] = mackskill(stability_measurement_envB,1);
    
    stability_measurement_sorted_by_days = {};
    stability_measurement_sorted_by_sessions = {};
    for mouse_ind = 1:number_of_mice
        stability_measurement_sorted_by_days{mouse_ind,1} = [stability_measurement_envA(mouse_ind,2:2:9);stability_measurement_envB(mouse_ind,1:4)]';
        stability_measurement_sorted_by_sessions{mouse_ind,1} = [stability_measurement_envA(mouse_ind,1:5);stability_measurement_envB(mouse_ind,:)]';
    end
    
    stability_measurement_sorted_by_days_all_mice = cell2mat(stability_measurement_sorted_by_days);
    stability_measurement_sorted_by_sessions_all_mice = cell2mat(stability_measurement_sorted_by_sessions);
    
    
    [~,time_stats] = mackskill(stability_measurement_sorted_by_days_all_mice,4);
    [~,experience_stats] = mackskill(stability_measurement_sorted_by_sessions_all_mice,5);
    
    time_exp_stats =  [drift_stats_envA.T,drift_stats_envA.df,drift_stats_envA.p;...
        drift_stats_envB.T,drift_stats_envB.df,drift_stats_envB.p;...
        time_stats.T,time_stats.df,time_stats.p;...
        experience_stats.T,experience_stats.df,experience_stats.p];
    
    
    col_names = {'Chi','df','pvalue'};
    row_names = {'Drift Env. A', 'Drift Env. B','Time','Experience'};
    statistics_table = array2table(time_exp_stats,'RowNames',row_names,'VariableNames',col_names);
    
elseif strcmp(flag,'context') % Related to Figures 3E-3G
    
    context_specificity_all_mice = stability_measurement;
    mean_context_specificity_across_mice = mean(context_specificity_all_mice,'omitnan');
    ste_context_specificity_across_mice = std(context_specificity_all_mice,'omitnan')./sqrt(number_of_mice);
    
    figure('units','normalized','position',[0.3 0.25 0.2 0.4])
    hold on
    plot(context_specificity_all_mice','color',[0.7 0.7 0.7])
    errorbar(mean_context_specificity_across_mice,ste_context_specificity_across_mice,'-o','color',[0.2 0.2 0.2],...
        'linewidth',2,'MarkerSize',4,'MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',[0.2 0.2 0.2])
    axis square
    ylabel(measurement_name)
    set(gca,'xtick',[1,2],'xticklabels',{'A-A','ABA'})
    xlim([0.5 2.5])
    ylim(ylims)
    title(figure_name,'FontSize',24,'Fontweight','normal')
    
    [p_exact,~,stats_exact] = signrank(context_specificity_all_mice(:,1),context_specificity_all_mice(:,2),'method','exact');
    [p_app,~,stats_app] = signrank(context_specificity_all_mice(:,1),context_specificity_all_mice(:,2),'method','approximate');
    
    context_stats = [stats_exact.signedrank,NaN,p_exact;...
        stats_app.signedrank,stats_app.zval,p_app];
    
    row_names = {'Exact','Approximated'};
    col_names = {'signed_ranks','Z','pvalue'};
    statistics_table = array2table(context_stats,'RowNames',row_names,'VariableNames',col_names);
    
elseif strcmp(flag,'controls') % Related to Figures S1A-S1F
    
    figure('units','normalized','position',[0.3 0.25 0.195 0.4])
    plt = [];
    for env_ind = 1:2
        current_days_ind = day_ind_per_env{env_ind};
        current_env_stability = cell2mat(stability_measurement(:,env_ind));
        
        mean_stability_across_mice = mean(current_env_stability,'omitnan');
        ste_stability_across_mice = std(current_env_stability,'omitnan')./sqrt(number_of_mice);
        
        hold on
        plt(env_ind) = errorbar(current_days_ind,mean_stability_across_mice,ste_stability_across_mice,'-o','color',env_colors(env_ind,:),'linewidth',2,'MarkerSize',4,...
            'MarkerEdgeColor',env_colors(env_ind,:),'MarkerFaceColor',env_colors(env_ind,:));
        
        if env_ind == 2
            legend(plt,{'Environment A','Environment B'},'Location','southwest')
            legend('boxoff')
        end
        
    end
    ylim(ylims)
    xlim([0 22])
    axis square
    ylabel(measurement_name)
    title(figure_name,'FontSize',24,'Fontweight','normal')
    xlabel('Day of experiment')
    axis square
    set(gca,'xtick',4:4:20)
    
    % sort data for statistical testing
    stability_measurement_envA = cell2mat(stability_measurement(:,1));
    stability_measurement_envB = cell2mat(stability_measurement(:,2));
    
    [~,drift_stats_envA] = mackskill(stability_measurement_envA,1);
    [~,drift_stats_envB] = mackskill(stability_measurement_envB,1);
    
    
    drift_stats =  [drift_stats_envA.T,drift_stats_envA.df,drift_stats_envA.p;...
        drift_stats_envB.T,drift_stats_envB.df,drift_stats_envB.p];
    
    
    col_names = {'Chi','df','pvalue'};
    row_names = {'Drift Env. A', 'Drift Env. B'};
    statistics_table = array2table(drift_stats,'RowNames',row_names,'VariableNames',col_names);
    
elseif strcmp(flag,'registration') % Related to Figure S1K
    
    elapsed_days_per_env_reg = {2:2:18,2:2:10,4:4:21;...
        2:2:18,4:4:16,4:4:21};
    
    elapsed_sessions_per_env_reg = {1:9,1:5,1:5;...
        1:9,2:2:8,1:5};
    
    env_colors_reg = [0.7, 0.7 0.7; env_colors];
    
    figure('units','normalized','position',[0.25 0.125 0.35 0.7])
    axes('pos',[0 1 1 1],'visible','off');
    text(0.425,-0.025,figure_name,'fontsize',24);
    
    time_exp_stats_both_measurements = {};
    for measurement_ind = 1:2
        
        current_measurement = stability_measurement{measurement_ind};
        current_measurement_name = measurement_name{measurement_ind};
        
        if measurement_ind == 1 % ensemble rate correlation
            plt_lgd = {'Environment A - all data','Environment A - first/last six sessions','Environment B'};
        elseif measurement_ind == 2 % tuning curve correlation
            plt_lgd = {'Environment A - all data','Environment A - odd/even days','Environment B'};
        end
        
        plt = [];
        for env_ind = 1:3
            elapsed_days = elapsed_days_per_env_reg{measurement_ind,env_ind};
            elapsed_sessions = elapsed_sessions_per_env_reg{measurement_ind,env_ind};
            current_env_stability = cell2mat(current_measurement(:,env_ind));
            
            number_of_mice = sum(~isnan(current_env_stability));
            
            mean_stability_across_mice = mean(current_env_stability,'omitnan');
            ste_stability_across_mice = std(current_env_stability,'omitnan')./sqrt(number_of_mice);
            
            subplot(2,2,1+2*(measurement_ind-1))
            hold on
            
            plt(env_ind) = errorbar(elapsed_days,mean_stability_across_mice, ste_stability_across_mice,'-o','color',env_colors_reg(env_ind,:),'linewidth',2,'MarkerSize',4,...
                'MarkerEdgeColor',env_colors_reg(env_ind,:),'MarkerFaceColor',env_colors_reg(env_ind,:));
            ylim(ylims)
            xlim([0 21])
            axis square
            ylabel(current_measurement_name)
            set(gca,'xtick',2:2:20)
            if env_ind == 3
                legend(plt,plt_lgd,'Location','southwest')
                legend('boxoff')
            end
            xlabel('Days difference')
            
            subplot(2,2,2+2*(measurement_ind-1))
            hold on
            errorbar(elapsed_sessions,mean_stability_across_mice, ste_stability_across_mice,'-o','color',env_colors_reg(env_ind,:),'linewidth',2,'MarkerSize',4,...
                'MarkerEdgeColor',env_colors_reg(env_ind,:),'MarkerFaceColor',env_colors_reg(env_ind,:));
            ylim(ylims)
            xlim([0 10])
            axis square
            xlabel('Sessions difference')
            set(gca,'xtick',1:9,'yticklabels',[])
        end
        
        
        % sort data for statistical testing
        stability_measurement_envA = cell2mat(current_measurement(:,2));
        stability_measurement_envB = cell2mat(current_measurement(:,3));
        
        
        stability_measurement_sorted_by_days = {};
        stability_measurement_sorted_by_sessions = {};
        for mouse_ind = 1:number_of_mice
            if measurement_ind == 1 % ensemble rate correlation
                stability_measurement_sorted_by_days{mouse_ind,1} = [stability_measurement_envA(mouse_ind,[2,4]);stability_measurement_envB(mouse_ind,1:2)]';
                stability_measurement_sorted_by_sessions{mouse_ind,1} = [stability_measurement_envA(mouse_ind,:);stability_measurement_envB(mouse_ind,1:5)]';
            elseif measurement_ind == 2 % tuning curve correlation
                stability_measurement_sorted_by_days{mouse_ind,1} = [stability_measurement_envA(mouse_ind,1:4);stability_measurement_envB(mouse_ind,1:4)]';
                stability_measurement_sorted_by_sessions{mouse_ind,1} = [stability_measurement_envA(mouse_ind,1:2);stability_measurement_envB(mouse_ind,[2,4])]';
            end
            
        end
        
        stability_measurement_sorted_by_days_all_mice = cell2mat(stability_measurement_sorted_by_days);
        stability_measurement_sorted_by_sessions_all_mice = cell2mat(stability_measurement_sorted_by_sessions);
        
        if measurement_ind == 1 % ensemble rate correlation
            [~,time_stats] = mackskill(stability_measurement_sorted_by_days_all_mice,2);
            [~,experience_stats] = mackskill(stability_measurement_sorted_by_sessions_all_mice,5);
        elseif measurement_ind == 2 % tuning curve correlation
            [~,time_stats] = mackskill(stability_measurement_sorted_by_days_all_mice,4);
            [~,experience_stats] = mackskill(stability_measurement_sorted_by_sessions_all_mice,2);
        end
        
        time_exp_stats =  [time_stats.T,time_stats.df,time_stats.p;...
            experience_stats.T,experience_stats.df,experience_stats.p];
        
        
        time_exp_stats_both_measurements{measurement_ind,1} = time_exp_stats;
        
    end
    
    col_names = {'Chi','df','pvalue'};
    row_names = {'Rate_Time','Rate_Experience','Tuning_Time','Tuning_Experience'};
    
    statistics_table = array2table(cell2mat(time_exp_stats_both_measurements),'RowNames',row_names,'VariableNames',col_names);
    
elseif strcmp(flag,'independence') % Related to Figure S1K
    
    load(['Y:\personal\nitzang\time_v_exp\all_2v4_days_mice\scripts\for_all_mice\magma_colormap.mat'])
    
    stability_measurement_pooled = cell2mat(stability_measurement(:));
    number_of_cells = size(stability_measurement_pooled,1);
    r = corr(stability_measurement_pooled(:,4),stability_measurement_pooled(:,ylims{1}));
    
    figure('units','normalized','position',[0.3 0.25 0.2 0.4])
    dscatter(stability_measurement_pooled(:,4),stability_measurement_pooled(:,ylims{1}))
    text(0.05,0.975,['R^2 = ',num2str(round(r^2,4))],'units','normalized')
    text(0.05,0.925,['N = ',num2str(number_of_cells)],'units','normalized')
    
    ylabel(measurement_name)
    title(figure_name,'FontSize',24,'Fontweight','normal')
    xlim([-1 1])
    ylim(ylims{2})
    colormap(magma_colormap)
    axis square
    ylabel('Absolute activity rate difference (Hz)')
    xlabel('Tuning curve correlation')
    
    cb = colorbar('SouthOutside');
    cb.Position = [0.6 0.8 0.3 0.03];
    cb.Ticks = [0.05 1];
    cb.YTickLabel = {'min','max'};


elseif strcmp(flag,'individual') % Related to Figure S2B-S2D
    
    axes('pos',[0 1 1 1],'visible','off');
    text(0.175,-0.035,figure_name,'fontsize',24);

    mouse_group_id = [1,1,1,2,2,2,1,1];
    environment_types = {'Env. A (Straight)','Env. B (L-shaped)'; 'Env. A (L-shaped)','Env. B (Straight)'};
    
    for mouse_ind = 1:number_of_mice
        current_mouse_stability_measurement = stability_measurement(mouse_ind,:);
        
        plt = [];
        for env_ind = 1:2
            elapsed_days = elapsed_days_per_env{env_ind};
            current_env_stability_measurement = current_mouse_stability_measurement{env_ind};
            
            subplot(number_of_mice,2,1+2*(mouse_ind-1))
            hold on
            plt(env_ind) = plot(elapsed_days, current_env_stability_measurement,'-','color',env_colors(env_ind,:),'linewidth',1.5);
            text(0.05,0.215-0.1*(env_ind-1),[environment_types{mouse_group_id(mouse_ind),env_ind}],...
                'units','normalized','color',env_colors(env_ind,:),'FontSize',6)
            
            ylim(ylims)
            xlim([0 21])
            ylabel(measurement_name,'FontSize',8)
            set(gca,'xtick',2:4:20,'ytick',ylims(1):0.2:ylims(2))
            
            if mouse_ind == 8
                xlabel('Days difference','FontSize',8)
            end
            
            
            subplot(number_of_mice,2,2+2*(mouse_ind-1))
            hold on
            plot( current_env_stability_measurement,'-','color',env_colors(env_ind,:),'linewidth',1.5);
            ylim(ylims)
            xlim([0 10])
            
            set(gca,'xtick',1:2:9,'yticklabels',[])
            
            if mouse_ind == 8
                xlabel('Sessions difference','FontSize',8)
            end
            
        end
        
    end
    
end

end
