function [u_n, v_n, time] = load_2DKSsolution(foldername, method, dt, T, N, L_s1, L_s2)

    mkdir([pwd  '/data/' foldername '' ]);

    phys_file = [pwd '/data/' foldername '/phys_' method '_N_' num2str(N) '' ...
        '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
    
    u_n = readmatrix(phys_file);

    
    four_file = [pwd '/data/' foldername '/four_' method '_N_' num2str(N) '' ...
        '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
    
    % v_n = readmatrix(four_file);
    v_n = u_n;

    time_file = [pwd '/data/' foldername '/time_' method '_N_' num2str(N) '' ...
        '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.dat'];
    
    time = readmatrix(time_file);

end