function save_measures(foldername, measure1, measure2, measure3, IC, N, dt, T, L_s1, L_s2)

    if not(isfolder([pwd  '/data' ]))                                           % create local directories for data storage
        mkdir([pwd  '/data/' foldername '_measures' ]);
    end
    
    switch foldername
        case 'energygrowth'
            L_scale = dt;
    
            l2mode_file = [pwd '/data/' foldername '_measures/L2mode_' strjoin(IC(1),'') '_compare' num2str(length(IC)) '_T_' ...
                num2str(T) '_start_' num2str(L_scale(1)) '_end_' num2str(L_scale(end)) '_target_' num2str(L_s1) '.dat'];
            
            writematrix(measure1, l2mode_file,'Delimiter','tab');
        
            l2avg_file = [pwd '/data/' foldername '_measures/L2avg_' strjoin(IC(1),'') '_compare' num2str(length(IC)) '_T_' ...
                num2str(T) '_start_' num2str(L_scale(1)) '_end_' num2str(L_scale(end)) '_target_' num2str(L_s1) '.dat'];
            
            writematrix(measure2, l2avg_file,'Delimiter','tab');
        
        case 'spatial'
        
            l2_file = [pwd '/data/' foldername '_measures/L2_' IC '_T_' ...
                num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
            
            writematrix(measure1, l2_file,'Delimiter','tab');
        
            inf_file = [pwd '/data/' foldername '_measures/LInf_' IC '_T_' ...
                num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
            
            writematrix(measure2, inf_file,'Delimiter','tab');
        
            comptime_file = [pwd '/data/' foldername '_measures/comptime_' IC '_T_' ...
                num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
            
            writematrix(measure3, comptime_file,'Delimiter','tab');
        case 'temporal'
        
            l2_file = [pwd '/data/' foldername '_measures/L2_' IC '_T_' ...
                num2str(T) '_N_' num2str(N) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
            
            writematrix(measure1, l2_file,'Delimiter','tab');
        
            inf_file = [pwd '/data/' foldername '_measures/LInf_' IC '_T_' ...
                num2str(T) '_N_' num2str(N) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
            
            writematrix(measure2, inf_file,'Delimiter','tab');
        
            comptime_file = [pwd '/data/' foldername '_measures/comptime_' IC '_T_' ...
                num2str(T) '_N_' num2str(N) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
            
            writematrix(measure3, comptime_file,'Delimiter','tab');
    
    end