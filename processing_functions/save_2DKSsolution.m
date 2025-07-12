function save_2DKSsolution(foldername, u_n, v_n, time, IC, dt, T, N, L_s1, L_s2, saved, originalIC)

    mkdir([pwd  '/data/' foldername '' ]);

    phys_file = [pwd '/data/' foldername '/phys_' IC '_N_' num2str(N) '' ...
    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples_' num2str(saved) '.dat'];
    four_file = [pwd '/data/' foldername '/four_' IC '_N_' num2str(N) '' ...
    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples_' num2str(saved) '.dat'];
    time_file = [pwd '/data/' foldername '/time_' IC '_N_' num2str(N) '' ...
    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples_' num2str(saved) '.dat'];
    switch IC
        case 'optimized'
            phys_file = [pwd '/data/' foldername '/phys_' IC '_' originalIC '_N_' num2str(N) '' ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples_' num2str(saved) '.dat'];
            four_file = [pwd '/data/' foldername '/four_' IC '_' originalIC '_N_' num2str(N) '' ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples_' num2str(saved) '.dat'];
            time_file = [pwd '/data/' foldername '/time_' IC '_' originalIC '_N_' num2str(N) '' ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples_' num2str(saved) '.dat'];
    end
    writematrix(u_n, phys_file,'Delimiter','tab');
    writematrix(v_n, four_file,'Delimiter','tab');
    writematrix(time, time_file,'Delimiter','tab');

end