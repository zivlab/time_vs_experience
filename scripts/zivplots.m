function [] = zivplots(tuning_curves,place_cells_ind,smoothing_sigma)

days_list = [2,4:4:20;1:4:21];
sessions_list = [1,2:2:10;1:6];
env_full_names = {'Environment A','Environment B'};

for env_ind = 1:2
    
    figure('units','normalized','position',[0.2+0.25*(env_ind-1) 0.15 0.2 0.55])
    axes('pos',[0 1 1 1],'visible','off');
    text(0.125,-0.035,env_full_names{env_ind},'fontsize',16);

    
    smoothed_activity_all_mice = {};
    pc_ind_all_mice = {};
    for mouse_ind = 1:8
        current_mouse_activity = tuning_curves{mouse_ind,env_ind};
        current_mouse_PC_ind = place_cells_ind{mouse_ind,env_ind};
        
        num_cells = size(current_mouse_activity,2)./2;
        for sess_ind = 1:size(current_mouse_activity,3)
         
            current_sess_activity = current_mouse_activity(:,:,sess_ind);
            pc_both_dir = current_mouse_PC_ind{sess_ind};
            
            current_sess_activity_smooth = zeros(size(current_sess_activity));
            for cell_id = 1:size(current_sess_activity,1)
                current_sess_activity_smooth(cell_id,:) = imgaussfilt(current_sess_activity(cell_id,:),smoothing_sigma);
            end
            smoothed_activity_all_mice{mouse_ind,sess_ind} = current_sess_activity_smooth;
            
            if mouse_ind > 1
                pc_ind_all_mice{mouse_ind,sess_ind} = pc_both_dir + size(cell2mat(smoothed_activity_all_mice(1:mouse_ind-1,1)),1);
            else
                pc_ind_all_mice{mouse_ind,sess_ind} = pc_both_dir;
            end
        end
        
    end
    
    num_sess = size(days_list,2);
    ind = 1;
    for sess1 = 1:num_sess
        
        ref_sess_activity = cell2mat(smoothed_activity_all_mice(:,sessions_list(env_ind,sess1)));
        ref_sess_pc_ind = cell2mat(pc_ind_all_mice(:,sessions_list(env_ind,sess1))');
        ref_sess_activity_pc =  ref_sess_activity(ref_sess_pc_ind,:);
        
        [sess_max,sess_peak] = max(ref_sess_activity_pc,[],2);
        [~,sess_sorted] = sort(sess_peak);
        ref_sess_activity_pc_norm = ref_sess_activity_pc./sess_max;
        
        for sess2 = 1:num_sess
            comp_sess_activity = cell2mat(smoothed_activity_all_mice(:,sessions_list(env_ind,sess2)));
            comp_sess_activity_pc = comp_sess_activity(ref_sess_pc_ind,:);
            sess_max = max(comp_sess_activity_pc,[],2);
            comp_sess_activity_pc_norm = comp_sess_activity_pc./sess_max;
            
            
            subplot(num_sess,num_sess,ind)
            imagesc(comp_sess_activity_pc_norm(sess_sorted,:),[0 1])
           
            if sess1 == 6
               xlabel(['Day ',num2str(days_list(env_ind,sess2))])
            end
            if sess2 == 1
               ylabel(['Day ',num2str(days_list(env_ind,sess1))])
            end
             
            if sess1 == 1 && sess2 == 6
                 title({'Position'},'units','normalized','Position',[0.5 1.01],'FontSize',10,'Fontweight','normal');  
            set(gca,'YAxisLocation','right')
             ylabel('Cell order','FontSize',10)
            end
              set(gca,'xtick', [],'ytick', []);
            ind = ind + 1;
            
        end
    end
    
    jet_colormap = colormap(jet);
    colormap(jet_colormap*0.9)
    
    cb = colorbar();
    cb.Position = [0.925 0.107 0.03 0.39];
    cb.Ticks = [0 1];
    ylabel(cb,'Normalized activity rate','FontSize',10,'Rotation',270);
end

end