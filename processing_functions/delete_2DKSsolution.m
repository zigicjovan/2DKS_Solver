function delete_2DKSsolution(foldername, IC, dt, T, N, L_s1, L_s2, Ntime_save_max)

    switch foldername
        case {'forward','backward'}

            phys_file = [pwd '/data/' foldername '/phys_' IC '_N_' num2str(N) '' ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples*'];
            
            delete(phys_file);

            four_file = [pwd '/data/' foldername '/four_' IC '_N_' num2str(N) '' ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples*'];
            
            delete(four_file);
        
            time_file = [pwd '/data/' foldername '/time_' IC '_N_' num2str(N) '' ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples*'];
            
            delete(time_file);

            try
                part = dt*Ntime_save_max;
                for Tpart = part:part:T

                        phys_file = [pwd '/data/' foldername '/phys_' IC '_N_' num2str(N) '' ...
                        '_T_' num2str(Tpart) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples*'];
                
                        delete(phys_file);
    
                        four_file = [pwd '/data/' foldername '/four_' IC '_N_' num2str(N) '' ...
                        '_T_' num2str(Tpart) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples*'];
                
                        delete(four_file);
            
                        time_file = [pwd '/data/' foldername '/time_' IC '_N_' num2str(N) '' ...
                        '_T_' num2str(Tpart) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_samples*'];
                
                        delete(time_file);
                end
            catch
            end
        case 'normL2'
            try
                norm_file = [pwd '/data/' foldername '/normL2_' IC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
                delete(norm_file);
            catch
                norm_file = [pwd '/data/' foldername '/normL2_' IC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '.dat'];
                delete(norm_file);
            end

        case 'normL2_t'
            norm_file = [pwd '/data/' foldername '/normL2_' IC '_N_' num2str(N) '' ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '.dat'];
            
            delete(norm_file);

        case 'spectrum'
            norm_file = [pwd '/data/' foldername '/spectrum_' IC '_N_' num2str(N) '' ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '.dat'];
            
            delete(norm_file);
        case {'optimization'}
            diagnostics_file = [pwd '/data/optimization/diagnostics_' IC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];

            linesearchJ_file = [pwd '/data/optimization/linesearchJ_' IC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];

            delete(diagnostics_file);
            delete(linesearchJ_file);
    end

end