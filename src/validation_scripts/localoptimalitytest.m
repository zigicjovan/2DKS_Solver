function hessapprox = localoptimalitytest(J_opt,u_IC_opt,IC,N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max)

% to test local maximum

    numsensitivitytests = 10;
    epsspace = (1./(2.^linspace(0,8,9)))';
    hessapprox = NaN(length(epsspace),numsensitivitytests);
    hessapproxnum = hessapprox;
    deltapert = hessapprox;
    cosangle = hessapprox;
    angle = hessapprox;
    u_opt_mag = sqrt(sum(abs(u_IC_opt(:)).^2)*(2*pi)^2*(L_s1*L_s2)/N^2);
    u_opt_sphere = sqrt(K)*(u_IC_opt/u_opt_mag);
    for j = 1:numsensitivitytests
        u_noise = (2*rand(N)-ones(N));
        u_pert_mag = sqrt(sum( u_noise(:) .* conj(u_noise(:)) )*(2*pi)^2*(L_s1*L_s2)/N^2);                  % compute norm of IC
        u_noise = u_noise/u_pert_mag;        
        IC_p = 'pertoptIC';
        for i = 1:length(epsspace)
            eps = epsspace(i,1);
            u_IC_optp = u_IC_opt + eps*u_noise(:); 
            [ v_TCp , ~ , ~ ] = solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,u_IC_optp,IC_p);
            pertmagplus = sum( abs(v_TCp(:)).^2 )*(2*pi)^2*(L_s1*L_s2)/(N*N)^2;
            pertdiffplus = pertmagplus - J_opt;
            u_IC_optp = u_IC_opt - eps*u_noise(:); 
            [ v_TCp , ~ , ~ ] = solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,u_IC_optp,IC_p);
            pertmagminus = sum( abs(v_TCp(:)).^2 )*(2*pi)^2*(L_s1*L_s2)/(N*N)^2;
            pertdiffminus = pertmagminus - J_opt;

            u_0p_mag = sqrt(sum( abs(u_IC_optp(:)).^2 )*(2*pi)^2*(L_s1*L_s2)/N^2);                  % compute norm of IC
            u_0p = sqrt(K)*(u_IC_optp/u_0p_mag); 
            deltapert(i,j) = sqrt(sum( abs(u_0p(:) - u_opt_sphere(:)).^2 )*(2*pi)^2*(L_s1*L_s2)/N^2);
            cosangle(i,j) = (sum(conj(u_0p(:)).*u_opt_sphere(:)) *(2*pi)^2*(L_s1*L_s2) / N^2) / K;   % should be cos(theta)
            angle(i,j) = acos(max(-1,min(1,real(cosangle(i,j)))));
            hessapproxnum(i,j) = pertdiffplus + pertdiffminus;
            hessapprox(i,j) = hessapproxnum(i,j)/(deltapert(i,j).^2);
        end
    end
end