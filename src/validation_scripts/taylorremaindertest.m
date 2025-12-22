function taylorremaindertest(J_opt,u_IC_opt,v_TC_opt,IC,N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,originalIC)

% ===============================
% Taylor remainder test for adjoint solver
% ===============================

numsensitivitytests = 5;
%epsspace = 1./(2.^linspace(0,8,9))';
epsspace = logspace(-14,-5,5)';

taylor_remainder = NaN(length(epsspace),numsensitivitytests);

% ---- gradient from adjoint ----
[~, gradJ, ~] = solve_2DKS( IC,'backward',N,K,L_s1,L_s2,dt,T, save_each,Ntime_save_max,v_TC_opt,originalIC);

for j = 1:numsensitivitytests

    % random perturbation
    v = 2*rand(N) - 1;
    v = v(:);

    % normalize direction
    vnorm = sqrt(sum(abs(v).^2)*(2*pi)^2*(L_s1*L_s2)/N^2);
    v = v / vnorm;

    % directional derivative
    dJ = real(sum(conj(gradJ(:)).*v(:))*(2*pi)^2*(L_s1*L_s2)/N^2);

    IC_p = 'pertoptIC';
    for i = 1:length(epsspace)
        eps = epsspace(i);

        % forward solve at perturbed IC
        u_p = u_IC_opt + eps*v;
        [v_TCp,~,~] = solve_2DKS( IC,'forward',N,K,L_s1,L_s2,dt,T, save_each,Ntime_save_max,u_p,IC_p);

        Jp = sum(abs(v_TCp(:)).^2)*(2*pi)^2*(L_s1*L_s2)/(N*N)^2;

        % Taylor remainder
        taylor_remainder(i,j) = ...
            abs(Jp - J_opt - eps*dJ) / eps^2;
    end
end

% ---- plot ----
figure; loglog(epsspace, taylor_remainder, 'o-');
xlabel('\epsilon'); ylabel('Taylor remainder');
title('Adjoint Taylor Test');
grid on;
