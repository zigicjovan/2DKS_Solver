%{
    mi = (i-1)*7;
    semilogy(l2norms_guess(:,mi+1),'r')
    hold on
    semilogy(l2norms_guess(:,mi+2),'g')
    semilogy(l2norms_guess(:,mi+3),'b')
    semilogy(l2norms_guess(:,mi+4),'m')
    semilogy(l2norms_opt(:,mi+1),'r')
    semilogy(l2norms_opt(:,mi+2),'g')
    semilogy(l2norms_opt(:,mi+3),'b')
    semilogy(l2norms_opt(:,mi+4),'m')


    for i = 1:4
    mi = (i-1)*7;
    semilogy(l2norms_guess(:,mi+0),'r--')
    hold on
    semilogy(l2norms_guess(:,mi+4),'r:')
    semilogy(l2norms_guess(:,mi+8),'b--')
    semilogy(l2norms_guess(:,mi+12),'b:')
    semilogy(l2norms_guess(:,mi+16),'g--')
    semilogy(l2norms_guess(:,mi+20),'g:')
    semilogy(l2norms_guess(:,mi+24),'m--')
    semilogy(l2norms_opt(:,mi+0),'r--')
    hold on
    semilogy(l2norms_opt(:,mi+4),'r:')
    semilogy(l2norms_opt(:,mi+8),'b--')
    semilogy(l2norms_opt(:,mi+12),'b:')
    semilogy(l2norms_opt(:,mi+16),'g--')
    semilogy(l2norms_opt(:,mi+20),'g:')
    semilogy(l2norms_opt(:,mi+24),'m--')
    %}

lamklist = NaN(numberoftests,1);
lamfitlist = lamklist;
residlist = lamklist;
counter = 0;
for lscale = 1:numberoftests
    %for guess = 1:4
        counter = counter + 1;
        ell = L_scale(lscale);
        pts = size(l2norms_opt,1); % /2 for half domain
        yl = log(l2norms_opt(1:end,counter)); % (1:length/2,counter) for half domain
        dx = linspace(0,T,pts); % T/2 for half domain
        xl = dx(1:end);
        J = [ ones(pts,1) , xl' ];
        x_fit = (J'*J)\J' * yl;
        y_fit = J*x_fit;
        resid_fit = sum((yl-y_fit).^2);
        if ell >= 1
            lamk = NaN(floor(ell)+1,floor(ell)+1);
            for k1 = 0:floor(ell)
                for k2 = 0:floor(ell)
                    lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
                end
            end
        else
            lamk = NaN(5,5);
            for k1 = 0:5
                for k2 = 0:5
                    lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
                end
            end
            lamk(1,1) = -inf;
        end
        lamklist(counter,1) = max(max(lamk));
        lamfitlist(counter,1) = x_fit(2,1);
        residlist(counter,1) = abs(lamfitlist(counter,1) - lamklist(counter,1));
    %end
end


time = linspace(0,T,length(l2norms_guess));
h = figure;
%{
semilogy(time,l2norms_opt(:,1),'LineWidth',0.5,'Marker','.')
hold on
semilogy(time,l2norms_opt(:,2),'LineWidth',0.5,'Marker','.')
semilogy(time,l2norms_opt(:,3),'LineWidth',0.5,'Marker','.')
%{
semilogy(time,l2norms_opt(:,4),'LineWidth',0.5,'Marker','.')
semilogy(time,l2norms_opt(:,5),'LineWidth',0.5,'Marker','.')
semilogy(time,l2norms_opt(:,6),'LineWidth',0.5,'Marker','.')
semilogy(time,l2norms_opt(:,7),'LineWidth',0.5,'Marker','.')
%}
%}
semilogy(time,l2norms_guess(:,1),'LineWidth',0.5,'Marker','.')
hold on
semilogy(time,l2norms_guess(:,2),'LineWidth',0.5,'Marker','.')
semilogy(time,l2norms_guess(:,3),'LineWidth',0.5,'Marker','.')
%{
semilogy(time,l2norms_guess(:,4),'LineWidth',0.5,'Marker','.')
semilogy(time,l2norms_guess(:,5),'LineWidth',0.5,'Marker','.')
semilogy(time,l2norms_guess(:,6),'LineWidth',0.5,'Marker','.')
semilogy(time,l2norms_guess(:,7),'LineWidth',0.5,'Marker','.')
%}
%}
set(gcf,'Position',[100 100 900 750])
xlabel('Time $t$','Interpreter','latex'); 
ylabel('$||\phi(t;\varphi)||_{L^2}$','Interpreter','latex');
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')    
%{
legend('$\ell = 1.1$','$\ell = 1.4$','$\ell = 1.5$','$\ell = 2.2$',...
    '$\ell = 3.2$','$\ell = 5.2$','$\ell = 10.2$',...
    'Interpreter','latex','Location','northwest')
%}
legend('$\ell = 0.5$','$\ell = 0.7$','$\ell = 0.9$', ...
    'Interpreter','latex','Location','southwest')
%}
%title('Evolution of optimized $L^2$ norm for varying isotropic domains','Interpreter','latex')
title('Evolution of $L^2$ norm for varying isotropic domains','Interpreter','latex')
subtitle(['$\varphi = \varphi_{' IC '}, T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
filename = [pwd '/media/optimization/validationguess_' IC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '.pdf'];
filenamefig = [pwd '/media/optimization/validationguess_' IC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '.fig'];
pause(0.5)
exportgraphics(h,filename)
saveas(h,filenamefig)