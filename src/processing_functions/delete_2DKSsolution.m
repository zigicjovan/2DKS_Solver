function delete_2DKSsolution(foldername, IC, dt, T, N, K, L_s1, L_s2, utility1, utility2)

    fullT = utility1(2);
    originalIC = utility2;
    parameterlist = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_T_' num2str(T) ];
    parameterlistanyT = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_t_*_T_' num2str(T) ];
    switch foldername
        case {'forward','backward'}
            folderpath = [pwd '/data/' foldername '/'];
            phys_file = [folderpath 'phys_' parameterlistanyT '_samples_*.bin'];
            four_file = [folderpath 'four_' parameterlistanyT '_samples_*.bin'];
            time_file = [folderpath 'time_' parameterlistanyT '_samples_*.bin'];
            switch IC
                case 'optimized'
                    phys_file = [folderpath 'phys_' originalIC '_' parameterlistanyT '_samples_*.bin'];
                    four_file = [folderpath 'four_' originalIC '_' parameterlistanyT '_samples_*.bin'];
                    time_file = [folderpath 'time_' originalIC '_' parameterlistanyT '_samples_*.bin'];
            end
            allphysfiles = dir(phys_file);
            allfourfiles = dir(four_file);
            alltimefiles = dir(time_file);
            numfiles = size(allphysfiles,1);
            for i = 1:numfiles
                cur_phys = fullfile(folderpath, allphysfiles(i).name);
                cur_four = fullfile(folderpath, allfourfiles(i).name);
                cur_time = fullfile(folderpath, alltimefiles(i).name);
                delete(cur_phys);   
                delete(cur_four);           
                delete(cur_time);
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