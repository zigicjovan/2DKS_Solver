function save_measures(foldername, measure1, measure2, measure3, IC, N, dt, T, K, L_s1, L_s2)

    switch foldername
        case {'energygrowth','spatial','temporal'}
            if not(isfolder([pwd  '/data/' foldername '_measures' ]))                                           % create local directories for data storage
                mkdir([pwd  '/data/' foldername '_measures' ]);
                addpath([pwd  '/data/' foldername '_measures' ]);
            end
    end

    parameterlist = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_T_' num2str(T) ];
    %originalIC = utility1;
    %tol = utility2;
    %optparameters = [ originalIC '_' parameterlist '_tol_' num2str(tol) ];

    switch foldername
        case 'optimization'
            Klist = unique(measure1(:,1));
            Llist = unique(measure1(:,2));
            Tlist = unique(measure1(:,3));
            for IC_i = 1:length(IC)
                IC_cur = strjoin(IC(IC_i));
                parameterlist = [IC_cur '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(Klist(1),'%.0f') '_' num2str(Klist(end),'%.0f') '_L_' num2str(Llist(1),'%.2f') '_' num2str(Llist(end),'%.2f') '_T_' num2str(Tlist(1),'%.2f') '_' num2str(Tlist(end),'%.2f') ];
                Jinit_file = [pwd '/data/' foldername '/Jinit_' parameterlist '.dat'];
                Jopt_file = [pwd '/data/' foldername '/Jopt_' parameterlist '.dat'];  
                writematrix([measure1(:,1:3),measure1(:,3+IC_i)], Jinit_file,'Delimiter','tab'); 
                writematrix([measure2(:,1:3),measure2(:,3+IC_i)], Jopt_file,'Delimiter','tab');
            end
        case 'energygrowth'
            L_scale = dt;
            l2mode_file = [pwd '/data/' foldername '_measures/L2mode_' strjoin(IC(1),'') '_compare' num2str(length(IC)) '_T_' ...
                num2str(T) '_start_' num2str(L_scale(1)) '_end_' num2str(L_scale(end)) '_target_' num2str(L_s1) '.dat'];
            l2avg_file = [pwd '/data/' foldername '_measures/L2avg_' strjoin(IC(1),'') '_compare' num2str(length(IC)) '_T_' ...
                num2str(T) '_start_' num2str(L_scale(1)) '_end_' num2str(L_scale(end)) '_target_' num2str(L_s1) '.dat'];
            writematrix(measure1, l2mode_file,'Delimiter','tab');            
            writematrix(measure2, l2avg_file,'Delimiter','tab');
        case 'spatial'
            l2_file = [pwd '/data/' foldername '_measures/L2_' parameterlist '.dat'];
            inf_file = [pwd '/data/' foldername '_measures/LInf_' parameterlist '.dat'];  
            comptime_file = [pwd '/data/' foldername '_measures/comptime_' parameterlist '.dat'];
            writematrix(measure1, l2_file,'Delimiter','tab'); 
            writematrix(measure2, inf_file,'Delimiter','tab');
            writematrix(measure3, comptime_file,'Delimiter','tab');
        case 'temporal'    
            l2_file = [pwd '/data/' foldername '_measures/L2_' parameterlist '.dat'];
            inf_file = [pwd '/data/' foldername '_measures/LInf_' parameterlist '.dat'];
            comptime_file = [pwd '/data/' foldername '_measures/comptime_' parameterlist '.dat'];
            writematrix(measure1, l2_file,'Delimiter','tab');      
            writematrix(measure2, inf_file,'Delimiter','tab');
            writematrix(measure3, comptime_file,'Delimiter','tab');
    end