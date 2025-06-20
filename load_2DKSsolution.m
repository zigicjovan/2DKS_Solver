function [file1, file2] = load_2DKSsolution(foldername, IC, dt, T, N, L_s1, L_s2)

    switch foldername
        case 'time_evolution'
            phys_file = [pwd '/data/' foldername '/phys_' IC '_N_' num2str(N) '' ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '.dat'];
            
            file1 = readmatrix(phys_file);
        
            time_file = [pwd '/data/' foldername '/time_' IC '_N_' num2str(N) '' ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '.dat'];
            
            file2 = readmatrix(time_file);
        case 'normL2'
            norm_file = [pwd '/data/' foldername '/normL2_' IC '_N_' num2str(N) '' ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '.dat'];
            
            file1 = readmatrix(norm_file);
            file2 = 0;
        case 'normL2_t'
            norm_file = [pwd '/data/' foldername '/normL2_' IC '_N_' num2str(N) '' ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '.dat'];
            
            file1 = readmatrix(norm_file);
            file2 = 0;
        case 'spectrum'
            norm_file = [pwd '/data/' foldername '/spectrum_' IC '_N_' num2str(N) '' ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '.dat'];
            
            file1 = readmatrix(norm_file);
            file2 = 0;
    end

end