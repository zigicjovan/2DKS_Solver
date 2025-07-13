function [file1, file2, file3] = load_2DKSsolution(foldername, IC, dt, T, N, L_s1, L_s2, saved, originalIC)

    switch foldername
        case {'forward','backward'}
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

            file1 = readmatrix(phys_file);   
            file2 = readmatrix(four_file);
            file3 = readmatrix(time_file);
        case 'normL2'
            try
                norm_file = [pwd '/data/' foldername '/normL2_' IC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
                switch IC
                    case 'optimized'
                        norm_file = [pwd '/data/' foldername '/normL2_' IC '_' originalIC '_N_' num2str(N) '' ...
                        '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
                end
                file1 = readmatrix(norm_file);
            catch
                norm_file = [pwd '/data/' foldername '/normL2_' IC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '.dat'];
                switch IC
                    case 'optimized'
                        norm_file = [pwd '/data/' foldername '/normL2_' IC '_' originalIC '_N_' num2str(N) '' ...
                        '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '.dat'];
                end
                file1 = readmatrix(norm_file);
            end
            
            file2 = 0;
            file3 = 0;
        case 'normL2_t'
            norm_file = [pwd '/data/' foldername '/normL2_' IC '_N_' num2str(N) '' ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '.dat'];
            switch IC
                case 'optimized'
                    norm_file = [pwd '/data/' foldername '/normL2_' IC '_' originalIC '_N_' num2str(N) '' ...
                    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '.dat'];
            end
            
            file1 = readmatrix(norm_file);
            file2 = 0;
            file3 = 0;
        case 'spectrum'
            norm_file = [pwd '/data/' foldername '/spectrum_' IC '_N_' num2str(N) '' ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '.dat'];
            switch IC
                case 'optimized'
                    norm_file = [pwd '/data/' foldername '/spectrum_' IC '_' originalIC '_N_' num2str(N) '' ...
                    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '.dat'];
            end
            
            file1 = readmatrix(norm_file);
            file2 = 0;
            file3 = 0;
        case {'optimization'}
            diagnostics_file = [pwd '/data/optimization/diagnostics_' IC '_' originalIC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];

            linesearchJ_file = [pwd '/data/optimization/linesearchJ_' IC '_' originalIC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];

            file1 = readmatrix(diagnostics_file);
            file2 = readmatrix(linesearchJ_file);
            file3 = 0;
    end

end