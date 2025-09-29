function delete_2DKSsolution(foldername, IC, dt, T, N, K, L_s1, L_s2, Ntime_save_max, originalIC)

    parameterlist = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_T_' num2str(T) ];
    switch foldername
        case {'forward','backward'}
            phys_file = [pwd '/data/' foldername '/phys_' parameterlist '_samples*'];
            four_file = [pwd '/data/' foldername '/four_' parameterlist '_samples*'];
            time_file = [pwd '/data/' foldername '/time_' parameterlist '_samples*'];
            switch IC
                case 'optimized'
                    phys_file = [pwd '/data/' foldername '/phys_' originalIC '_' parameterlist '_samples*'];
                    four_file = [pwd '/data/' foldername '/four_' originalIC '_' parameterlist '_samples*'];
                    time_file = [pwd '/data/' foldername '/time_' originalIC '_' parameterlist '_samples*'];
            end
            delete(phys_file);   
            delete(four_file);           
            delete(time_file);
            try
                part = dt*Ntime_save_max;
                for Tpart = part:part:T
                        phys_file = [pwd '/data/' foldername '/phys_' parameterlist '_samples*'];
                        four_file = [pwd '/data/' foldername '/four_' parameterlist '_samples*'];
                        time_file = [pwd '/data/' foldername '/time_' parameterlist '_samples*'];
                        switch IC
                            case 'optimized'
                                phys_file = [pwd '/data/' foldername '/phys_' originalIC '_' parameterlist '_samples*'];
                                four_file = [pwd '/data/' foldername '/four_' originalIC '_' parameterlist '_samples*'];
                                time_file = [pwd '/data/' foldername '/time_' originalIC '_' parameterlist '_samples*'];
                        end
                        delete(phys_file);
                        delete(four_file);
                        delete(time_file);
                end
            catch
            end
        case 'normL2'
            try
                norm_file = [pwd '/data/' foldername '/normL2_' parameterlist '.dat'];
                delete(norm_file);
            catch
                norm_file = [pwd '/data/' foldername '/normL2_' parameterlist '.dat'];
                switch IC
                    case 'optimized'
                        norm_file = [pwd '/data/' foldername '/normL2_' originalIC '_' parameterlist '.dat'];
                end
                delete(norm_file);
            end

        case 'normL2_t'
            norm_file = [pwd '/data/' foldername '/normL2_' parameterlist '.dat'];
            switch IC
                case 'optimized'
                    norm_file = [pwd '/data/' foldername '/normL2_' originalIC '_' parameterlist '.dat'];
            end
            delete(norm_file);

        case 'spectrum'
            norm_file = [pwd '/data/' foldername '/spectrum_' parameterlist '.dat'];
            switch IC
                case 'optimized'
                    norm_file = [pwd '/data/' foldername '/spectrum_' originalIC '_' parameterlist '.dat'];
            end
            delete(norm_file);
        case {'optimization'}
            diagnostics_file = [pwd '/data/optimization/diagnostics_' originalIC '_' parameterlist '.dat'];
            linesearchJ_file = [pwd '/data/optimization/linesearchJ_' originalIC '_' parameterlist '.dat'];
            delete(diagnostics_file);
            delete(linesearchJ_file);
    end

end