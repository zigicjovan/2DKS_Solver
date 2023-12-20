function save_2DKSsolution(foldername, u_n, v_n, time, method, dt, T, N, L_s1, L_s2)

    mkdir([pwd  '/data/' foldername '' ]);

    phys_file = [pwd '/data/' foldername '/phys_' method '_N_' num2str(N) '' ...
        '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
    
    writematrix(u_n, phys_file,'Delimiter','tab');

    four_file = [pwd '/data/' foldername '/four_' method '_N_' num2str(N) '' ...
        '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
    
    writematrix(v_n, four_file,'Delimiter','tab');

    time_file = [pwd '/data/' foldername '/time_' method '_N_' num2str(N) '' ...
        '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
    
    writematrix(time, time_file,'Delimiter','tab');

end