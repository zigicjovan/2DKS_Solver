function [file1, file2, file3] = load_2DKSsolution(foldername, IC, dt, T, N, K, L_s1, L_s2, utility1, utility2)

    parameterlist = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_T_' num2str(T) ];
    saved = utility1;
    originalIC = utility2;
    switch foldername
        case {'forward','backward'}
            phys_file = [pwd '/data/' foldername '/phys_' parameterlist '_samples_' num2str(saved) '.dat'];
            four_file = [pwd '/data/' foldername '/four_' parameterlist '_samples_' num2str(saved) '.dat'];
            time_file = [pwd '/data/' foldername '/time_' parameterlist '_samples_' num2str(saved) '.dat'];
            switch IC
                case 'optimized'
                    phys_file = [pwd '/data/' foldername '/phys_' originalIC '_' parameterlist '_samples_' num2str(saved) '.dat'];
                    four_file = [pwd '/data/' foldername '/four_' originalIC '_' parameterlist '_samples_' num2str(saved) '.dat'];
                    time_file = [pwd '/data/' foldername '/time_' originalIC '_' parameterlist '_samples_' num2str(saved) '.dat'];
            end
            file1 = readmatrix(phys_file);   
            file2 = readmatrix(four_file);
            file3 = readmatrix(time_file);
        case {'optimal'}
            tol = utility1;
            phys_file = [pwd '/data/' foldername '/physIC_' parameterlist '_tol_' num2str(tol) '.dat'];
            four_file = [pwd '/data/' foldername '/fourTC_' parameterlist '_tol_' num2str(tol) '.dat'];
            time_file = [pwd '/data/' foldername '/time_' parameterlist '_tol_' num2str(tol) '.dat'];
            file1 = readmatrix(phys_file);   
            file2 = readmatrix(four_file);
            file3 = readmatrix(time_file);
        case 'energyL2'
            tol = utility1;
            norm_file = [pwd '/data/' foldername '/energyL2_' parameterlist '.dat'];
            switch IC
                case 'optimized'
                    norm_file = [pwd '/data/' foldername '/energyL2_' originalIC '_' parameterlist '_tol_' num2str(tol) '.dat'];
            end
            file1 = readmatrix(norm_file);
            file2 = 0;
            file3 = 0;
        case 'energyL2_t'
            tol = utility1;
            norm_file = [pwd '/data/' foldername '/energyL2_' parameterlist '.dat'];
            switch IC
                case 'optimized'
                    norm_file = [pwd '/data/' foldername '/energyL2_' originalIC '_' parameterlist '_tol_' num2str(tol) '.dat'];
            end         
            file1 = readmatrix(norm_file);
            file2 = 0;
            file3 = 0;
        case 'spectrum'
            tol = utility1;
            norm_file = [pwd '/data/' foldername '/spectrum_' parameterlist '.dat'];
            switch IC
                case 'optimized'
                    norm_file = [pwd '/data/' foldername '/spectrum_' originalIC '_' parameterlist '_tol_' num2str(tol) '.dat'];
            end            
            file1 = readmatrix(norm_file);
            file2 = 0;
            file3 = 0;
        case {'optimization'}
            tol = utility1;
            diagnostics_file = [pwd '/data/optimization/diagnostics_' originalIC '_' parameterlist '_tol_' num2str(tol) '.dat'];
            linesearchJ_file = [pwd '/data/optimization/linesearchJ_' originalIC '_' parameterlist '_tol_' num2str(tol) '.dat'];
            file1 = readmatrix(diagnostics_file);
            file2 = readmatrix(linesearchJ_file);
            file3 = 0;
    end

end