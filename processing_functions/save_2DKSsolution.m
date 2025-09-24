function save_2DKSsolution(foldername, u_n, v_n, time, IC, dt, T, N, L_s1, L_s2, saved, utility1)

    if not(isfolder([pwd  '/data/' foldername '' ]))                                           % create local directories for data storage
        mkdir([pwd  '/data/' foldername '' ]);
        addpath([pwd  '/data/' foldername '' ]);
    end

    if strcmp('optimal',foldername)
        tol = utility1;
        phys_file = [pwd '/data/' foldername '/phys_' IC '_N_' num2str(N) '' ...
        '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_tol_' num2str(tol) '.dat'];
        four_file = [pwd '/data/' foldername '/four_' IC '_N_' num2str(N) '' ...
        '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_tol_' num2str(tol) '.dat'];
        time_file = [pwd '/data/' foldername '/time_' IC '_N_' num2str(N) '' ...
        '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_tol_' num2str(tol) '.dat'];
        writematrix(u_n, phys_file);
        writematrix(v_n, four_file);
        writematrix(time, time_file);
    else
        phys_file = [pwd '/data/' foldername '/phys_' IC '_N_' num2str(N) '' ...
        '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples_' num2str(saved) '.dat'];
        four_file = [pwd '/data/' foldername '/four_' IC '_N_' num2str(N) '' ...
        '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples_' num2str(saved) '.dat'];
        time_file = [pwd '/data/' foldername '/time_' IC '_N_' num2str(N) '' ...
        '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples_' num2str(saved) '.dat'];
        switch IC
            case 'optimized'
                originalIC = utility1;
                phys_file = [pwd '/data/' foldername '/phys_' IC '_' originalIC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples_' num2str(saved) '.dat'];
                four_file = [pwd '/data/' foldername '/four_' IC '_' originalIC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples_' num2str(saved) '.dat'];
                time_file = [pwd '/data/' foldername '/time_' IC '_' originalIC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples_' num2str(saved) '.dat'];
        end
        writematrix(u_n, phys_file);
        writematrix(v_n, four_file);
        writematrix(time, time_file);
    end
end