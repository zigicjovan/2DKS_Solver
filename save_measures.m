function save_measures(foldername, error_2, error_inf, comptime, method, N, dt, T, L_s1, L_s2)

switch foldername
    case 'spatial'
        mkdir([pwd  '/data/' foldername '_measures' ]);
    
        l2_file = [pwd '/data/' foldername '_measures/L2_' method '_T_' ...
            num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
        
        writematrix(error_2, l2_file,'Delimiter','tab');
    
        inf_file = [pwd '/data/' foldername '_measures/LInf_' method '_T_' ...
            num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
        
        writematrix(error_inf, inf_file,'Delimiter','tab');
    
        comptime_file = [pwd '/data/' foldername '_measures/comptime_' method '_T_' ...
            num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
        
        writematrix(comptime, comptime_file,'Delimiter','tab');
    case 'temporal'
        mkdir([pwd  '/data/' foldername '_measures' ]);
    
        l2_file = [pwd '/data/' foldername '_measures/L2_' method '_T_' ...
            num2str(T) '_N_' num2str(N) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
        
        writematrix(error_2, l2_file,'Delimiter','tab');
    
        inf_file = [pwd '/data/' foldername '_measures/LInf_' method '_T_' ...
            num2str(T) '_N_' num2str(N) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
        
        writematrix(error_inf, inf_file,'Delimiter','tab');
    
        comptime_file = [pwd '/data/' foldername '_measures/comptime_' method '_T_' ...
            num2str(T) '_N_' num2str(N) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
        
        writematrix(comptime, comptime_file,'Delimiter','tab');

end