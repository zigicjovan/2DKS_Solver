function delete_2DKSsolution(foldername, IC, dt, T, N, K, L_s1, L_s2, utility1, utility2)

    Ntime_save_max = utility1(1);
    fullT = utility1(2);
    originalIC = utility2;
    parameterlistT = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_t_' num2str(T) '_T_' num2str(fullT) ];
    parameterlist = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_T_' num2str(T) ];
    switch foldername
        case {'forward','backward'}
            phys_file = [pwd '/data/' foldername '/phys_' parameterlistT '_samples*'];
            four_file = [pwd '/data/' foldername '/four_' parameterlistT '_samples*'];
            time_file = [pwd '/data/' foldername '/time_' parameterlistT '_samples*'];
            switch IC
                case 'optimized'
                    phys_file = [pwd '/data/' foldername '/phys_' originalIC '_' parameterlistT '_samples*'];
                    four_file = [pwd '/data/' foldername '/four_' originalIC '_' parameterlistT '_samples*'];
                    time_file = [pwd '/data/' foldername '/time_' originalIC '_' parameterlistT '_samples*'];
            end
            delete(phys_file);   
            delete(four_file);           
            delete(time_file);
            try
                part = dt*Ntime_save_max;
                for Tpart = part:part:T
                        parameterlistT = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_t_' num2str(Tpart) '_T_' num2str(fullT) ];
                        phys_file = [pwd '/data/' foldername '/phys_' parameterlistT '_samples*'];
                        four_file = [pwd '/data/' foldername '/four_' parameterlistT '_samples*'];
                        time_file = [pwd '/data/' foldername '/time_' parameterlistT '_samples*'];
                        switch IC
                            case 'optimized'
                                phys_file = [pwd '/data/' foldername '/phys_' originalIC '_' parameterlistT '_samples*'];
                                four_file = [pwd '/data/' foldername '/four_' originalIC '_' parameterlistT '_samples*'];
                                time_file = [pwd '/data/' foldername '/time_' originalIC '_' parameterlistT '_samples*'];
                        end
                        delete(phys_file);
                        delete(four_file);
                        delete(time_file);
                end
            catch
            end
        case 'energyL2'
            try
                norm_file = [pwd '/data/' foldername '/energyL2_' parameterlist '.dat'];
                delete(norm_file);
            catch
                norm_file = [pwd '/data/' foldername '/energyL2_' parameterlist '.dat'];
                switch IC
                    case 'optimized'
                        norm_file = [pwd '/data/' foldername '/energyL2_' originalIC '_' parameterlist '.dat'];
                end
                delete(norm_file);
            end

        case {'energyL2_t','spectrum'}
            norm_file = [pwd '/data/' foldername '/' foldername '_' parameterlist '.dat'];
            switch IC
                case 'optimized'
                    norm_file = [pwd '/data/' foldername '/' foldername '_' originalIC '_' parameterlist '.dat'];
            end
            delete(norm_file);

        case {'optimization'}
            diagnostics_file = [pwd '/data/optimization/diagnostics_' originalIC '_' parameterlist '.dat'];
            linesearchJ_file = [pwd '/data/optimization/linesearchJ_' originalIC '_' parameterlist '.dat'];
            delete(diagnostics_file);
            delete(linesearchJ_file);
    end

end