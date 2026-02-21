function [file1, file2, file3] = load_2DKSsolution(foldername, IC, dt, T, N, K, L_s1, L_s2, utility1, utility2)

    fullT = utility1(2);
    parameterlistT = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_t_' num2str(T) '_T_' num2str(fullT) ];
    parameterlist = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_T_' num2str(T) ];
    parameterlistcont = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_T_' ];
    saved = utility1(1);
    originalIC = utility2;
    switch foldername
        case {'forward','backward'}
            %phys_file = [pwd '/data/' foldername '/phys_' parameterlistT '_samples_' num2str(saved) '.bin'];
            four_file = [pwd '/data/' foldername '/four_' parameterlistT '_samples_' num2str(saved) '.bin'];
            %time_file = [pwd '/data/' foldername '/time_' parameterlistT '_samples_' num2str(saved) '.bin'];
            switch IC
                case 'optimized'
                    %phys_file = [pwd '/data/' foldername '/phys_' originalIC '_' parameterlistT '_samples_' num2str(saved) '.bin'];
                    four_file = [pwd '/data/' foldername '/four_' originalIC '_' parameterlistT '_samples_' num2str(saved) '.bin'];
                    %time_file = [pwd '/data/' foldername '/time_' originalIC '_' parameterlistT '_samples_' num2str(saved) '.bin'];
            end
            %file1 = read_binary(phys_file,N,N,false);   
            file1 = 0;
            file2 = read_binary(four_file,N,N,true);
            %file3 = read_binary(time_file,1,1,false);
            file3 = 0;
        case {'optimal'}
            tol = utility1(1);
            % use exact or nearest previous T
            try
                phys_file = sprintf('%s/data/%s/physIC_%s_tol_%g.bin', pwd, foldername, parameterlist, tol);
                four_file = sprintf('%s/data/%s/fourTC_%s_tol_%g.bin', pwd, foldername, parameterlist, tol);
                time_file = sprintf('%s/data/%s/time_%s_tol_%g.bin', pwd, foldername, parameterlist, tol);
                file1 = read_binary(phys_file,N,N,false);   
                file2 = read_binary(four_file,N,N,true);
                file3 = read_binary(time_file,1,1,false);
            catch
                phys_filegroup = sprintf('%s/data/%s/physIC_%s*_tol_%g.bin', pwd, foldername, parameterlistcont, tol);
                phys_files = dir(phys_filegroup);
                
                numfiles = numel(phys_files);
                Tvals = zeros(numfiles,1);
                for k = 1:numfiles
                    tokens = regexp(phys_files(k).name, ...
                        '_T_([0-9eE\+\-\.]+)_tol_', ...
                        'tokens', 'once');
                
                    Tvals(k) = str2double(tokens{1});
                end
                mask = Tvals < (T+1e-10);
                Tvals = Tvals(mask);
                prevT = Tvals(end);
                phys_file = sprintf('%s/data/%s/physIC_%s%g_tol_%g.bin', pwd, foldername, parameterlistcont, prevT, tol);
                four_file = sprintf('%s/data/%s/fourTC_%s%g_tol_%g.bin', pwd, foldername, parameterlistcont, prevT, tol);
                time_file = sprintf('%s/data/%s/time_%s%g_tol_%g.bin', pwd, foldername, parameterlistcont, prevT, tol);
                file1 = read_binary(phys_file,N,N,false);   
                file2 = read_binary(four_file,N,N,true);
                file3 = read_binary(time_file,1,1,false);
            end
        case {'energyL2','energyL2_t','spectrum'}
            tol = utility1(1);
            norm_file = [pwd '/data/' foldername '/' foldername '_' parameterlist '.dat'];
            if strcmp(IC, 'optimized')
                norm_file = [pwd '/data/' foldername '/' foldername '_' originalIC '_' parameterlist '_tol_' num2str(tol) '.dat'];
            end
            file1 = readmatrix(norm_file);
            file2 = 0;
            file3 = 0;
        case {'optimization'}
            tol = utility1(1);
            diagnostics_file = [pwd '/data/optimization/diagnostics/diagnostics_' originalIC '_' parameterlist '_tol_' num2str(tol) '.dat'];
            linesearchJ_file = [pwd '/data/optimization/linesearch/linesearchJ_' originalIC '_' parameterlist '_tol_' num2str(tol) '.dat'];
            file1 = readmatrix(diagnostics_file);
            file2 = readmatrix(linesearchJ_file);
            file3 = 0;
    end

end