function save_2DKSsolution(foldername, u_n, v_n, time, IC, dt, T, N, K, L_s1, L_s2, utility1, utility2)

    if not(isfolder([pwd  '/data/' foldername '' ]))                                           % create local directories for data storage
        mkdir([pwd  '/data/' foldername '' ]);
        addpath([pwd  '/data/' foldername '' ]);
    end

    fullT = utility1(2);
    saved = utility1(1);
    parameterlistT = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_t_' num2str(T) '_T_' num2str(fullT) ];
    parameterlist = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_T_' num2str(T) ];
    
    if strcmp('optimal',foldername)
        tol = utility2;
        phys_file = [pwd '/data/' foldername '/physIC_' parameterlist '_tol_' num2str(tol) '.dat'];
        four_file = [pwd '/data/' foldername '/fourTC_' parameterlist '_tol_' num2str(tol) '.dat'];
        time_file = [pwd '/data/' foldername '/time_' parameterlist '_tol_' num2str(tol) '.dat'];
        writematrix(u_n, phys_file);
        writematrix(v_n, four_file);
        writematrix(time, time_file);
    else
        phys_file = [pwd '/data/' foldername '/phys_' parameterlistT '_samples_' num2str(saved) '.dat'];
        four_file = [pwd '/data/' foldername '/four_' parameterlistT '_samples_' num2str(saved) '.dat'];
        time_file = [pwd '/data/' foldername '/time_' parameterlistT '_samples_' num2str(saved) '.dat'];
        switch IC
            case 'optimized'
                originalIC = utility2;
                phys_file = [pwd '/data/' foldername '/phys_' originalIC '_' parameterlistT '_samples_' num2str(saved) '.dat'];
                four_file = [pwd '/data/' foldername '/four_' originalIC '_' parameterlistT '_samples_' num2str(saved) '.dat'];
                time_file = [pwd '/data/' foldername '/time_' originalIC '_' parameterlistT '_samples_' num2str(saved) '.dat'];
        end
        writematrix(u_n, phys_file);
        writematrix(v_n, four_file);
        writematrix(time, time_file);
    end
end