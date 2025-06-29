function save_2DKSsolution(foldername, u_n, v_n, time, IC, dt, T, N, L_s1, L_s2, saved)

    mkdir([pwd  '/data/' foldername '' ]);

    phys_file = [pwd '/data/' foldername '/phys_' IC '_N_' num2str(N) '' ...
    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples_' num2str(saved) '.dat'];
    
    writematrix(u_n, phys_file,'Delimiter','tab');

    four_file = [pwd '/data/' foldername '/four_' IC '_N_' num2str(N) '' ...
    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples_' num2str(saved) '.dat'];
            
    writematrix(v_n, four_file,'Delimiter','tab');

    time_file = [pwd '/data/' foldername '/time_' IC '_N_' num2str(N) '' ...
    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples_' num2str(saved) '.dat'];
    
    writematrix(time, time_file,'Delimiter','tab');

end