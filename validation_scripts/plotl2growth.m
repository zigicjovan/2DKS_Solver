%% Test 1: macheps linear growth F %%

T = 150;
N = 512;
dt = 1e-3;

% 1.1

u_n = u11fn14;
ell1=1.1;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu11 = log(normu);
dx = linspace(0,T,pts);
normu11 = normu;

% validate growth rate
samplestart = 50;
sampleend = 100;
xl = dx(1,samplestart:sampleend);

%%% least squares fit %%%
yl = lognormu11( samplestart:sampleend , 1 );
J = [ ones(length(yl),1) , xl' ];
x_fit = (J'*J)\J' * yl;
y_fit = J*x_fit;
resid_fit = sum((yl-y_fit).^2);
ell = ell1;
lamk = NaN(floor(ell)+1,floor(ell)+1);
for k1 = 0:floor(ell)
    for k2 = 0:floor(ell)
        lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
    end
end
lam11 = max(max(lamk));
resid11 = norm(x_fit(2,1) - lam11) / norm(lam11);
xl11 = xl;
lamk11 = lamk;
lamfit11 = x_fit(2,1);

% 1.4

u_n = u14fn14;
ell1=1.4;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu14 = log(normu);
dx = linspace(0,T,pts);
normu14 = normu;

% validate growth rate
samplestart = 50;
sampleend = 100;
xl = dx(1,samplestart:sampleend);

%%% least squares fit %%%
yl = lognormu14( samplestart:sampleend , 1 );
J = [ ones(length(yl),1) , xl' ];
x_fit = (J'*J)\J' * yl;
y_fit = J*x_fit;
resid_fit = sum((yl-y_fit).^2);
ell = ell1;
lamk = NaN(floor(ell)+1,floor(ell)+1);
for k1 = 0:floor(ell)
    for k2 = 0:floor(ell)
        lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
    end
end
lam14 = max(max(lamk));
resid14 = norm(x_fit(2,1) - lam14) / norm(lam14);
xl14 = xl;
lamk14 = lamk;
lamfit14 = x_fit(2,1);

% 1.5

u_n = u15fn14;
ell1=1.5;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu15 = log(normu);
dx = linspace(0,T,pts);
normu15 = normu;

% validate growth rate
samplestart = 50;
sampleend = 100;
xl = dx(1,samplestart:sampleend);

%%% least squares fit %%%
yl = lognormu15( samplestart:sampleend , 1 );
J = [ ones(length(yl),1) , xl' ];
x_fit = (J'*J)\J' * yl;
y_fit = J*x_fit;
resid_fit = sum((yl-y_fit).^2);
ell = ell1;
lamk = NaN(floor(ell)+1,floor(ell)+1);
for k1 = 0:floor(ell)
    for k2 = 0:floor(ell)
        lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
    end
end
lam15 = max(max(lamk));
resid15 = norm(x_fit(2,1) - lam15) / norm(lam15);
xl15 = xl;
lamk15 = lamk;
lamfit15 = x_fit(2,1);

% 2.2

u_n = u22fn14;
ell1=2.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu22 = log(normu);
dx = linspace(0,T,pts);
normu22 = normu;

% validate growth rate
samplestart = 50;
sampleend = 100;
xl = dx(1,samplestart:sampleend);

%%% least squares fit %%%
yl = lognormu22( samplestart:sampleend , 1 );
J = [ ones(length(yl),1) , xl' ];
x_fit = (J'*J)\J' * yl;
y_fit = J*x_fit;
resid_fit = sum((yl-y_fit).^2);
ell = ell1;
lamk = NaN(floor(ell)+1,floor(ell)+1);
for k1 = 0:floor(ell)
    for k2 = 0:floor(ell)
        lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
    end
end
lam22 = max(max(lamk));
resid22 = norm(x_fit(2,1) - lam22) / norm(lam22);
xl22 = xl;
lamk22 = lamk;
lamfit22 = x_fit(2,1);

% 3.2

u_n = u32fn14;
ell1=3.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu32 = log(normu);
dx = linspace(0,T,pts);
normu32 = normu;

% validate growth rate
samplestart = 50;
sampleend = 100;
xl = dx(1,samplestart:sampleend);

%%% least squares fit %%%
yl = lognormu32( samplestart:sampleend , 1 );
J = [ ones(length(yl),1) , xl' ];
x_fit = (J'*J)\J' * yl;
y_fit = J*x_fit;
resid_fit = sum((yl-y_fit).^2);
ell = ell1;
lamk = NaN(floor(ell)+1,floor(ell)+1);
for k1 = 0:floor(ell)
    for k2 = 0:floor(ell)
        lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
    end
end
lam32 = max(max(lamk));
resid32 = norm(x_fit(2,1) - lam32) / norm(lam32);
xl32 = xl;
lamk32 = lamk;
lamfit32 = x_fit(2,1);

% 5.2

u_n = u52fn14;
ell1=5.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu52 = log(normu);
dx = linspace(0,T,pts);
normu52 = normu;

% validate growth rate
samplestart = 50;
sampleend = 100;
xl = dx(1,samplestart:sampleend);

%%% least squares fit %%%
yl = lognormu52( samplestart:sampleend , 1 );
J = [ ones(length(yl),1) , xl' ];
x_fit = (J'*J)\J' * yl;
y_fit = J*x_fit;
resid_fit = sum((yl-y_fit).^2);
ell = ell1;
lamk = NaN(floor(ell)+1,floor(ell)+1);
for k1 = 0:floor(ell)
    for k2 = 0:floor(ell)
        lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
    end
end
lam52 = max(max(lamk));
resid52 = norm(x_fit(2,1) - lam52) / norm(lam52);
xl52 = xl;
lamk52 = lamk;
lamfit52 = x_fit(2,1);

% 10.2

u_n = u102fn14;
ell1=10.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu102 = log(normu);
dx = linspace(0,T,pts);
normu102 = normu;

% validate growth rate
samplestart = 50;
sampleend = 100;
xl = dx(1,samplestart:sampleend);

%%% least squares fit %%%
yl = lognormu102( samplestart:sampleend , 1 );
J = [ ones(length(yl),1) , xl' ];
x_fit = (J'*J)\J' * yl;
y_fit = J*x_fit;
resid_fit = sum((yl-y_fit).^2);
ell = ell1;
lamk = NaN(floor(ell)+1,floor(ell)+1);
for k1 = 0:floor(ell)
    for k2 = 0:floor(ell)
        lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
    end
end
lam102 = max(max(lamk));
resid102 = norm(x_fit(2,1) - lam102) / norm(lam102);
xl102 = xl;
lamk102 = lamk;
lamfit102 = x_fit(2,1);

%figure(2)
semilogy(dx,normu11,'r-+','LineWidth',1)
hold on;
semilogy(dx,normu14,'g-+','LineWidth',1)
semilogy(dx,normu15,'b-+','LineWidth',1)
semilogy(dx,normu22,'r-o','LineWidth',1)
semilogy(dx,normu32,'g-o','LineWidth',1)
semilogy(dx,normu52,'b-o','LineWidth',1)
semilogy(dx,normu102,'r-x','LineWidth',1)
semilogy(xl11,10^(-16.9)*exp(lam11.*xl11),'LineWidth',3)
semilogy(xl14,10^(-16.6)*exp(lam14.*xl14),'LineWidth',3)
semilogy(xl15,10^(-16.75)*exp(lam15.*xl15),'LineWidth',3)
semilogy(xl22,10^(-16.9)*exp(lam22.*xl22),'LineWidth',3)
semilogy(xl32,10^(-16.45)*exp(lam32.*xl32),'LineWidth',3)
semilogy(xl52,10^(-16.4)*exp(lam52.*xl52),'LineWidth',3)
semilogy(xl102,10^(-16.0)*exp(lam102.*xl102),'LineWidth',3)

set(gcf,'Position',[100 100 900 750])
xlabel('$t$','Interpreter','latex'); 
ylabel('$|| \phi(t) ||_{L^2}$','Interpreter','latex');
xlim([0 T])
ylim([1e-17 1e0])

fontsize(12,"points")
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')    
title("2DKS growth rate for varying isotropic domains")
legend("$\ell = 1.1$", "$\ell = 1.4$",  "$\ell = 1.5$",  "$\ell = 2.2$",  "$\ell = 3.2$",  ...
    "$\ell = 5.2$",  "$\ell = 10.2$","$\lambda^*(\ell = 1.1)$", "$\lambda^*(\ell = 1.4)$", "$\lambda^*(\ell = 1.5)$", ...
    "$\lambda^*(\ell = 2.2)$", "$\lambda^*(\ell = 3.2)$", "$\lambda^*(\ell = 5.2)$",  "$\lambda^*(\ell = 10.2)$", ...
                'Location','northwest','NumColumns',2,'Interpreter','latex')

%% Test 2: sinL linear decay F %%

T = 150;
N = 512;
dt = 1e-3;

% 0.5

u_n = u05fn14;
ell1=0.5;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu05 = log(normu);
dx = linspace(0,T,pts);
normu05 = normu;

% validate growth rate
samplestart = 1;
sampleend = 3;
xl = dx(1,samplestart:sampleend);

%%% least squares fit %%%
yl = lognormu05( samplestart:sampleend , 1 );
J = [ ones(length(yl),1) , xl' ];
x_fit = (J'*J)\J' * yl;
y_fit = J*x_fit;
resid_fit = sum((yl-y_fit).^2);
ell = ell1;
lamk = NaN(5,5);
for k1 = 0:5
    for k2 = 0:5
        lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
    end
end
lamk(1,1) = -inf;
lam05 = max(max(lamk));
resid05 = norm(x_fit(2,1) - lam05) / norm(lam05);
xl05 = xl;
lamk05 = lamk;
lamfit05 = x_fit(2,1);

% 0.7

u_n = u07fn14;
ell1=0.7;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu07 = log(normu);
dx = linspace(0,T,pts);
normu07 = normu;

% validate growth rate
samplestart = 5;
sampleend = 15;
xl = dx(1,samplestart:sampleend);

%%% least squares fit %%%
yl = lognormu07( samplestart:sampleend , 1 );
J = [ ones(length(yl),1) , xl' ];
x_fit = (J'*J)\J' * yl;
y_fit = J*x_fit;
resid_fit = sum((yl-y_fit).^2);
ell = ell1;
lamk = NaN(5,5);
for k1 = 0:5
    for k2 = 0:5
        lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
    end
end
lamk(1,1) = -inf;
lam07 = max(max(lamk));
resid07 = norm(x_fit(2,1) - lam07) / norm(lam07);
xl07 = xl;
lamk07 = lamk;
lamfit07 = x_fit(2,1);

% 0.9

u_n = u09fn14;
ell1=0.9;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu09 = log(normu);
dx = linspace(0,T,pts);
normu09 = normu;

% validate growth rate
samplestart = 5;
sampleend = 15;
xl = dx(1,samplestart:sampleend);

%%% least squares fit %%%
yl = lognormu09( samplestart:sampleend , 1 );
J = [ ones(length(yl),1) , xl' ];
x_fit = (J'*J)\J' * yl;
y_fit = J*x_fit;
resid_fit = sum((yl-y_fit).^2);
ell = ell1;
lamk = NaN(5,5);
for k1 = 0:5
    for k2 = 0:5
        lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
    end
end
lamk(1,1) = -inf;
lam09 = max(max(lamk));
resid09 = norm(x_fit(2,1) - lam09) / norm(lam09);
xl09 = xl;
lamk09 = lamk;
lamfit09 = x_fit(2,1);

%%% Decay %%%
semilogy(dx(1,1:20),normu05(1:20,1),'-+','LineWidth',2)
hold on;
semilogy(dx(1,1:20),normu07(1:20,1),'-+','LineWidth',2)
semilogy(dx(1,1:20),normu09(1:20,1),'-+','LineWidth',2)
semilogy(xl05,10^(0)*exp(lam05.*xl05),'-','LineWidth',3)
semilogy(xl07,10^(0)*exp(lam07.*xl07),'-','LineWidth',3)
semilogy(xl09,10^(0)*exp(lam09.*xl09),'-','LineWidth',3)

set(gcf,'Position',[100 100 900 750])
xlabel('$t$','Interpreter','latex'); 
ylabel('$|| \phi(t) ||_{L^2}$','Interpreter','latex');
xlim([0 15])
ylim([1e-17 2e0])

fontsize(12,"points")
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')    
title("2DKS decay rate for varying isotropic domains")
legend("$\ell = 0.5$", "$\ell = 0.7$",  "$\ell = 0.9$", ...
    "$\lambda^*(\ell = 0.5)$",  "$\lambda^*(\ell = 0.7)$", "$\lambda^*(\ell = 0.9)$",  ...
                'Location','southeast','NumColumns',2,'Interpreter','latex')

%% Test 3: mn1,5,10 nonlinear growth MF %%

T = 150;
N = 512;
dt = 1e-3;

%%% fortran data %%%

% 1.1

u_n = u11fn1;
ell1=1.1;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu11fn1 = normu;

u_n = u11fn5;
ell1=1.1;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu11fn5 = normu;

u_n = u11fn10;
ell1=1.1;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu11fn10 = normu;

% 1.4

u_n = u14fn1;
ell1=1.4;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu14fn1 = normu;

u_n = u14fn5;
ell1=1.4;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu14fn5 = normu;

u_n = u14fn10;
ell1=1.4;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu14fn10 = normu;


% 1.5

u_n = u15fn1;
ell1=1.5;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu15fn1 = normu;

u_n = u15fn5;
ell1=1.5;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu15fn5 = normu;

u_n = u15fn10;
ell1=1.5;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu15fn10 = normu;

% 2.2

u_n = u22fn1;
ell1=2.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu22fn1 = normu;

u_n = u22fn5;
ell1=2.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu22fn5 = normu;

u_n = u22fn10;
ell1=2.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu22fn10 = normu;

% 3.2

u_n = u32fn1;
ell1=3.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu32fn1 = normu;

u_n = u32fn5;
ell1=3.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu32fn5 = normu;

u_n = u32fn10;
ell1=3.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu32fn10 = normu;

% 5.2

u_n = u52fn1;
ell1=5.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu52fn1 = normu;

u_n = u52fn5;
ell1=5.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu52fn5 = normu;

u_n = u52fn10;
ell1=5.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu52fn10 = normu;

% 10.2

u_n = u102fn1;
ell1=10.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu102fn1 = normu;

u_n = u102fn5;
ell1=10.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu102fn5 = normu;

u_n = u102fn10;
ell1=10.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu102fn10 = normu;

%%% Load Matlab data %%%

T = 150;
N = 512;
dt = 1e-3;
xx = 1.1;
yy = xx;
[u11n1, ~, ~] = load_2DKSsolution('time_evolution', 'mn1', 1e-3, T, N, xx, yy);
[u11n5, ~, ~] = load_2DKSsolution('time_evolution', 'mn5', 1e-3, T, N, xx, yy);
[u11n10, ~, ~] = load_2DKSsolution('time_evolution', 'mn10', 1e-3, T, N, xx, yy);
xx = 1.4;
yy = xx;
[u14n1, ~, ~] = load_2DKSsolution('time_evolution', 'mn1', 1e-3, T, N, xx, yy);
[u14n5, ~, ~] = load_2DKSsolution('time_evolution', 'mn5', 1e-3, T, N, xx, yy);
[u14n10, ~, ~] = load_2DKSsolution('time_evolution', 'mn10', 1e-3, T, N, xx, yy);
xx = 1.5;
yy = xx;
[u15n1, ~, ~] = load_2DKSsolution('time_evolution', 'mn1', 1e-3, T, N, xx, yy);
[u15n5, ~, ~] = load_2DKSsolution('time_evolution', 'mn5', 1e-3, T, N, xx, yy);
[u15n10, ~, ~] = load_2DKSsolution('time_evolution', 'mn10', 1e-3, T, N, xx, yy);
xx = 2.2;
yy = xx;
[u22n1, ~, ~] = load_2DKSsolution('time_evolution', 'mn1', 1e-3, T, N, xx, yy);
[u22n5, ~, ~] = load_2DKSsolution('time_evolution', 'mn5', 1e-3, T, N, xx, yy);
[u22n10, ~, ~] = load_2DKSsolution('time_evolution', 'mn10', 1e-3, T, N, xx, yy);
xx = 3.2;
yy = xx;
[u32n1, ~, ~] = load_2DKSsolution('time_evolution', 'mn1', 1e-3, T, N, xx, yy);
[u32n5, ~, ~] = load_2DKSsolution('time_evolution', 'mn5', 1e-3, T, N, xx, yy);
[u32n10, ~, ~] = load_2DKSsolution('time_evolution', 'mn10', 1e-3, T, N, xx, yy);
xx = 5.2;
yy = xx;
[u52n1, ~, ~] = load_2DKSsolution('time_evolution', 'mn1', 1e-3, T, N, xx, yy);
[u52n5, ~, ~] = load_2DKSsolution('time_evolution', 'mn5', 1e-3, T, N, xx, yy);
[u52n10, ~, ~] = load_2DKSsolution('time_evolution', 'mn10', 1e-3, T, N, xx, yy);
xx = 10.2;
yy = xx;
[u102n1, ~, ~] = load_2DKSsolution('time_evolution', 'mn1', 1e-3, T, N, xx, yy);
[u102n5, ~, ~] = load_2DKSsolution('time_evolution', 'mn5', 1e-3, T, N, xx, yy);
[u102n10, ~, ~] = load_2DKSsolution('time_evolution', 'mn10', 1e-3, T, N, xx, yy);


%% Matlab data %%%


% 1.1

u_n = u11n1;
ell1=1.1;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu11n1 = normu;

u_n = u11n5;
ell1=1.1;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu11n5 = normu;

u_n = u11n10;
ell1=1.1;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu11n10 = normu;

% 1.4

u_n = u14n1;
ell1=1.4;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu14n1 = normu;

u_n = u14n5;
ell1=1.4;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu14n5 = normu;

u_n = u14n10;
ell1=1.4;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu14n10 = normu;


% 1.5

u_n = u15n1;
ell1=1.5;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu15n1 = normu;

u_n = u15n5;
ell1=1.5;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu15n5 = normu;

u_n = u15n10;
ell1=1.5;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu15n10 = normu;

% 2.2

u_n = u22n1;
ell1=2.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu22n1 = normu;

u_n = u22n5;
ell1=2.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu22n5 = normu;

u_n = u22n10;
ell1=2.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu22n10 = normu;

% 3.2

u_n = u32n1;
ell1=3.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu32n1 = normu;

u_n = u32n5;
ell1=3.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu32n5 = normu;

u_n = u32n10;
ell1=3.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu32n10 = normu;

% 5.2

u_n = u52n1;
ell1=5.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu52n1 = normu;

u_n = u52n5;
ell1=5.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu52n5 = normu;

u_n = u52n10;
ell1=5.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu52n10 = normu;

% 10.2

u_n = u102n1;
ell1=10.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu102n1 = normu;

u_n = u102n5;
ell1=10.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu102n5 = normu;

u_n = u102n10;
ell1=10.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu102n10 = normu;

% + o x for 1 5 10 , r2 g2 b2 k2 r1 g1 b1 for 11 14 15 22 32 52 102 
figure(1)
semilogy(dx(1, 1:101),normu11fn1(1:101, 1),'r-+','LineWidth',2)
hold on;
semilogy(dx(1, 1:101),normu14fn1(1:101, 1),'g-+','LineWidth',2)
hold on;
semilogy(dx(1, 1:101),normu15fn1(1:101, 1),'b-+','LineWidth',2)
semilogy(dx(1, 1:101),normu22fn1(1:101, 1),'k-+','LineWidth',2)
semilogy(dx(1, 1:101),normu32fn1(1:101, 1),'r-+','LineWidth',1)
semilogy(dx(1, 1:101),normu52fn1(1:101, 1),'g-+','LineWidth',1)
semilogy(dx(1, 1:101),normu102fn1(1:101, 1),'b-+','LineWidth',1)

semilogy(dx(1, 1:101),normu11fn5(1:101, 1),'r-o','LineWidth',2)
semilogy(dx(1, 1:101),normu14fn5(1:101, 1),'g-o','LineWidth',2)
semilogy(dx(1, 1:101),normu15fn5(1:101, 1),'b-o','LineWidth',2)
semilogy(dx(1, 1:101),normu22fn5(1:101, 1),'k-o','LineWidth',2)
semilogy(dx(1, 1:101),normu32fn5(1:101, 1),'r-o','LineWidth',1)
semilogy(dx(1, 1:101),normu52fn5(1:101, 1),'g-o','LineWidth',1)
semilogy(dx(1, 1:101),normu102fn5(1:101, 1),'b-o','LineWidth',1)

semilogy(dx(1, 1:101),normu11fn10(1:101, 1),'r-x','LineWidth',2)
semilogy(dx(1, 1:101),normu14fn10(1:101, 1),'g-x','LineWidth',2)
semilogy(dx(1, 1:101),normu15fn10(1:101, 1),'b-x','LineWidth',2)
semilogy(dx(1, 1:101),normu22fn10(1:101, 1),'k-x','LineWidth',2)
semilogy(dx(1, 1:101),normu32fn10(1:101, 1),'r-x','LineWidth',1)
semilogy(dx(1, 1:101),normu52fn10(1:101, 1),'g-x','LineWidth',1)
semilogy(dx(1, 1:101),normu102fn10(1:101, 1),'b-x','LineWidth',1)

semilogy(dx(1, 1:101),normu11n1(1:101, 1),'r-+','LineWidth',2)
hold on;
semilogy(dx(1, 1:101),normu14n1(1:101, 1),'g-+','LineWidth',2)
semilogy(dx(1, 1:101),normu15n1(1:101, 1),'b-+','LineWidth',2)
semilogy(dx(1, 1:101),normu22n1(1:101, 1),'k-+','LineWidth',2)
semilogy(dx(1, 1:101),normu32n1(1:101, 1),'r-+','LineWidth',1)
semilogy(dx(1, 1:101),normu52n1(1:101, 1),'g-+','LineWidth',1)
semilogy(dx(1, 1:101),normu102n1(1:101, 1),'b-+','LineWidth',1)

semilogy(dx(1, 1:101),normu11n5(1:101, 1),'r-o','LineWidth',2)
semilogy(dx(1, 1:101),normu14n5(1:101, 1),'g-o','LineWidth',2)
semilogy(dx(1, 1:101),normu15n5(1:101, 1),'b-o','LineWidth',2)
semilogy(dx(1, 1:101),normu22n5(1:101, 1),'k-o','LineWidth',2)
semilogy(dx(1, 1:101),normu32n5(1:101, 1),'r-o','LineWidth',1)
semilogy(dx(1, 1:101),normu52n5(1:101, 1),'g-o','LineWidth',1)
semilogy(dx(1, 1:101),normu102n5(1:101, 1),'b-o','LineWidth',1)

semilogy(dx(1, 1:101),normu11n10(1:101, 1),'r-x','LineWidth',2)
semilogy(dx(1, 1:101),normu14n10(1:101, 1),'g-x','LineWidth',2)
semilogy(dx(1, 1:101),normu15n10(1:101, 1),'b-x','LineWidth',2)
semilogy(dx(1, 1:101),normu22n10(1:101, 1),'k-x','LineWidth',2)
semilogy(dx(1, 1:101),normu32n10(1:101, 1),'r-x','LineWidth',1)
semilogy(dx(1, 1:101),normu52n10(1:101, 1),'g-x','LineWidth',1)
semilogy(dx(1, 1:101),normu102n10(1:101, 1),'b-x','LineWidth',1)
%}
set(gcf,'Position',[100 100 900 750])
xlabel('$t$','Interpreter','latex'); 
ylabel('$|| \phi(t) ||_{L^2}$','Interpreter','latex');

fontsize(12,"points")
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')    
title("2DKS energy evolution for varying isotropic domains")
legend( "$\ell = 1.4$, IC n1", "$\ell = 1.5$, IC n1", ...
    "$\ell = 2.2$, IC n1","$\ell = 3.2$, IC n1", "$\ell = 5.2$, IC n1", "$\ell = 10.2$, IC n1",...
    "$\ell = 1.4$, IC n5", "$\ell = 1.5$, IC n5", ...
    "$\ell = 2.2$, IC n5","$\ell = 3.2$, IC n5", "$\ell = 5.2$, IC n5", "$\ell = 10.2$, IC n5",...
    "$\ell = 1.4$, IC n10", "$\ell = 1.5$, IC n10", ...
    "$\ell = 2.2$, IC n10","$\ell = 3.2$, IC n10", "$\ell = 5.2$, IC n10", "$\ell = 10.2$, IC n10",...
                'Location','southoutside','NumColumns',3,'Interpreter','latex')

%% relative difference (experimental: MATLAB; theoretical: fortran)
diffL2_11n1 = NaN(pts,1);
diffL2_14n1 = NaN(pts,1);
diffL2_15n1 = NaN(pts,1);
diffL2_22n1 = NaN(pts,1);
diffL2_32n1 = NaN(pts,1);
diffL2_52n1 = NaN(pts,1);
diffL2_102n1 = NaN(pts,1);
diffL2_11n5 = NaN(pts,1);
diffL2_14n5 = NaN(pts,1);
diffL2_15n5 = NaN(pts,1);
diffL2_22n5 = NaN(pts,1);
diffL2_32n5 = NaN(pts,1);
diffL2_52n5 = NaN(pts,1);
diffL2_102n5 = NaN(pts,1);
diffL2_11n10 = NaN(pts,1);
diffL2_14n10 = NaN(pts,1);
diffL2_15n10 = NaN(pts,1);
diffL2_22n10 = NaN(pts,1);
diffL2_32n10 = NaN(pts,1);
diffL2_52n10 = NaN(pts,1);
diffL2_102n10 = NaN(pts,1);
for j = 1:pts
    diffL2_11n1(j,1) = norm(normu11n1(j,1)-normu11fn1(j,1))/normu11fn1(j,1);
    diffL2_14n1(j,1) = norm(normu14n1(j,1)-normu14fn1(j,1))/normu14fn1(j,1);
    diffL2_15n1(j,1) = norm(normu15n1(j,1)-normu15fn1(j,1))/normu15fn1(j,1);
    diffL2_22n1(j,1) = norm(normu22n1(j,1)-normu22fn1(j,1))/normu22fn1(j,1);
    diffL2_32n1(j,1) = norm(normu32n1(j,1)-normu32fn1(j,1))/normu32fn1(j,1);
    diffL2_52n1(j,1) = norm(normu52n1(j,1)-normu52fn1(j,1))/normu52fn1(j,1);
    diffL2_102n1(j,1) = norm(normu102n1(j,1)-normu102fn1(j,1))/normu102fn1(j,1);
    
    diffL2_11n5(j,1) = norm(normu11n5(j,1)-normu11fn5(j,1))/normu11fn5(j,1);
    diffL2_14n5(j,1) = norm(normu14n5(j,1)-normu14fn5(j,1))/normu14fn5(j,1);
    diffL2_15n5(j,1) = norm(normu15n5(j,1)-normu15fn5(j,1))/normu15fn5(j,1);
    diffL2_22n5(j,1) = norm(normu22n5(j,1)-normu22fn5(j,1))/normu22fn5(j,1);
    diffL2_32n5(j,1) = norm(normu32n5(j,1)-normu32fn5(j,1))/normu32fn5(j,1);
    diffL2_52n5(j,1) = norm(normu52n5(j,1)-normu52fn5(j,1))/normu52fn5(j,1);
    diffL2_102n5(j,1) = norm(normu102n5(j,1)-normu102fn5(j,1))/normu102fn5(j,1);
        
    diffL2_11n10(j,1) = norm(normu11n10(j,1)-normu11fn10(j,1))/normu11fn10(j,1);
    diffL2_14n10(j,1) = norm(normu14n10(j,1)-normu14fn10(j,1))/normu14fn10(j,1);
    diffL2_15n10(j,1) = norm(normu15n10(j,1)-normu15fn10(j,1))/normu15fn10(j,1);
    diffL2_22n10(j,1) = norm(normu22n10(j,1)-normu22fn10(j,1))/normu22fn10(j,1);
    diffL2_32n10(j,1) = norm(normu32n10(j,1)-normu32fn10(j,1))/normu32fn10(j,1);
    diffL2_52n10(j,1) = norm(normu52n10(j,1)-normu52fn10(j,1))/normu52fn10(j,1);
    diffL2_102n10(j,1) = norm(normu102n10(j,1)-normu102fn10(j,1))/normu102fn10(j,1);
end
for j = 1
    if diffL2_11n1(j,1) == 0
        diffL2_11n1(j,1) = 1e-15;
    end
    if diffL2_14n1(j,1) == 0
        diffL2_14n1(j,1) = 1e-15;
    end
    if diffL2_15n1(j,1) == 0
        diffL2_15n1(j,1) = 1e-15;
    end
    if diffL2_22n1(j,1) == 0
        diffL2_22n1(j,1) = 1e-15;
    end
    if diffL2_32n1(j,1) == 0
        diffL2_32n1(j,1) = 1e-15;
    end
    if diffL2_52n1(j,1) == 0
        diffL2_52n1(j,1) = 1e-15;
    end
    if diffL2_102n1(j,1) == 0
        diffL2_102n1(j,1) = 1e-15;
    end
    
    if diffL2_11n5(j,1) == 0
        diffL2_11n5(j,1) = 1e-15;
    end
    if diffL2_14n5(j,1) == 0
        diffL2_14n5(j,1) = 1e-15;
    end
    if diffL2_15n5(j,1) == 0
        diffL2_15n5(j,1) = 1e-15;
    end
    if diffL2_22n5(j,1) == 0
        diffL2_22n5(j,1) = 1e-15;
    end
    if diffL2_32n5(j,1) == 0
        diffL2_32n5(j,1) = 1e-15;
    end
    if diffL2_52n5(j,1) == 0
        diffL2_52n5(j,1) = 1e-15;
    end
    if diffL2_102n5(j,1) == 0
        diffL2_102n5(j,1) = 1e-15;
    end

    if diffL2_11n10(j,1) == 0
        diffL2_11n10(j,1) = 1e-15;
    end
    if diffL2_14n10(j,1) == 0
        diffL2_14n10(j,1) = 1e-15;
    end
    if diffL2_15n10(j,1) == 0
        diffL2_15n10(j,1) = 1e-15;
    end
    if diffL2_22n10(j,1) == 0
        diffL2_22n10(j,1) = 1e-15;
    end
    if diffL2_32n10(j,1) == 0
        diffL2_32n10(j,1) = 1e-15;
    end
    if diffL2_52n10(j,1) == 0
        diffL2_52n10(j,1) = 1e-15;
    end
    if diffL2_102n10(j,1) == 0
        diffL2_102n10(j,1) = 1e-15;
    end
end

% + o x for 1 5 10 , r2 g2 b2 k2 r1 g1 b1 for 11 14 15 22 32 52 102 
figure(2)
semilogy(dx(1, 1:101),diffL2_11n1(1:101, 1),'r-+','LineWidth',2)

semilogy(dx(1, 1:101),diffL2_14n1(1:101, 1),'g-+','LineWidth',2)
hold on;
semilogy(dx(1, 1:101),diffL2_15n1(1:101, 1),'b-+','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_22n1(1:101, 1),'k-+','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_32n1(1:101, 1),'r-+','LineWidth',1)
semilogy(dx(1, 1:101),diffL2_52n1(1:101, 1),'g-+','LineWidth',1)
semilogy(dx(1, 1:101),diffL2_102n1(1:101, 1),'b-+','LineWidth',1)

semilogy(dx(1, 1:101),diffL2_11n5(1:101, 1),'r-o','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_14n5(1:101, 1),'g-o','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_15n5(1:101, 1),'b-o','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_22n5(1:101, 1),'k-o','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_32n5(1:101, 1),'r-o','LineWidth',1)
semilogy(dx(1, 1:101),diffL2_52n5(1:101, 1),'g-o','LineWidth',1)
semilogy(dx(1, 1:101),diffL2_102n5(1:101, 1),'b-o','LineWidth',1)

semilogy(dx(1, 1:101),diffL2_11n10(1:101, 1),'r-x','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_14n10(1:101, 1),'g-x','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_15n10(1:101, 1),'b-x','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_22n10(1:101, 1),'k-x','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_32n10(1:101, 1),'r-x','LineWidth',1)
semilogy(dx(1, 1:101),diffL2_52n10(1:101, 1),'g-x','LineWidth',1)
semilogy(dx(1, 1:101),diffL2_102n10(1:101, 1),'b-x','LineWidth',1)

set(gcf,'Position',[100 100 900 750])
        xlabel('$t$','Interpreter','latex'); 
        ylabel('$\frac{|| \phi_F(t) - \phi_M(t) ||_{L^2} }{ || \phi_F(t) ||_{L^2}}$','Interpreter','latex');

        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("Relative difference, MATLAB v. Fortran")
        legend( "$\ell = 1.4$, IC n1", "$\ell = 1.5$, IC n1", ...
    "$\ell = 2.2$, IC n1","$\ell = 3.2$, IC n1", "$\ell = 5.2$, IC n1", "$\ell = 10.2$, IC n1",...
    "$\ell = 1.4$, IC n5", "$\ell = 1.5$, IC n5", ...
    "$\ell = 2.2$, IC n5","$\ell = 3.2$, IC n5", "$\ell = 5.2$, IC n5", "$\ell = 10.2$, IC n5",...
    "$\ell = 1.4$, IC n10", "$\ell = 1.5$, IC n10", ...
    "$\ell = 2.2$, IC n10","$\ell = 3.2$, IC n10", "$\ell = 5.2$, IC n10", "$\ell = 10.2$, IC n10",...
                'Location','southoutside','NumColumns',3,'Interpreter','latex')
%}


%% Test 4: sinL MFK %%

% load matlab data
T = 100;
N = 128;
dt = 1e-3;
[u11sinLmat, ~, ~] = load_2DKSsolution('time_evolution', 'sinL', 1e-3, T, N, 1.1, 1.1);
[u12sinLmat, ~, ~] = load_2DKSsolution('time_evolution', 'sinL', 1e-3, T, N, 1.2, 1.2);
[u13sinLmat, ~, ~] = load_2DKSsolution('time_evolution', 'sinL', 1e-3, T, N, 1.3, 1.3);
[u14sinLmat, ~, ~] = load_2DKSsolution('time_evolution', 'sinL', 1e-3, T, N, 1.4, 1.4);
[u15sinLmat, ~, ~] = load_2DKSsolution('time_evolution', 'sinL', 1e-3, T, N, 1.5, 1.5);
[u16sinLmat, ~, ~] = load_2DKSsolution('time_evolution', 'sinL', 1e-3, T, N, 1.6, 1.6);
[u17sinLmat, ~, ~] = load_2DKSsolution('time_evolution', 'sinL', 1e-3, T, N, 1.7, 1.7);
[u18sinLmat, ~, ~] = load_2DKSsolution('time_evolution', 'sinL', 1e-3, T, N, 1.8, 1.8);
[u19sinLmat, ~, ~] = load_2DKSsolution('time_evolution', 'sinL', 1e-3, T, N, 1.9, 1.9);

% fortran

N = 512;
u_n = u11sinL;
ell1=1.1;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu11sinL = normu;

u_n = u12sinL;
ell1=1.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu12sinL = normu;

u_n = u13sinL;
ell1=1.3;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu13sinL = normu;

u_n = u14sinL;
ell1=1.4;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu14sinL = normu;

u_n = u15sinL;
ell1=1.5;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu15sinL = normu;

u_n = u16sinL;
ell1=1.6;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu16sinL = normu;

u_n = u17sinL;
ell1=1.7;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu17sinL = normu;

u_n = u18sinL;
ell1=1.8;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu18sinL = normu;

u_n = u19sinL;
ell1=1.9;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu19sinL = normu;

% matlab

N = 128;
u_n = u11sinLmat;
ell1=1.1;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu11sinLmat = normu;

u_n = u12sinLmat;
ell1=1.2;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu12sinLmat = normu;

u_n = u13sinLmat;
ell1=1.3;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu13sinLmat = normu;

u_n = u14sinLmat;
ell1=1.4;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu14sinLmat = normu;

u_n = u15sinLmat;
ell1=1.5;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu15sinLmat = normu;

u_n = u16sinLmat;
ell1=1.6;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu16sinLmat = normu;

u_n = u17sinLmat;
ell1=1.7;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu17sinLmat = normu;

u_n = u18sinLmat;
ell1=1.8;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu18sinLmat = normu;

u_n = u19sinLmat;
ell1=1.9;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu19sinLmat = normu;

u_filename = dir('udata_X1.1Y1.1_N256.*').name;
u11k = readmatrix(u_filename());
u_filename = dir('udata_X1.2Y1.2_N256.*').name;
u12k = readmatrix(u_filename());
u_filename = dir('udata_X1.3Y1.3_N256.*').name;
u13k = readmatrix(u_filename());
u_filename = dir('udata_X1.4Y1.4_N256.*').name;
u14k = readmatrix(u_filename());
u_filename = dir('udata_X1.5Y1.5_N256.*').name;
u15k = readmatrix(u_filename());
u_filename = dir('udata_X1.6Y1.6_N256.*').name;
u16k = readmatrix(u_filename());
u_filename = dir('udata_X1.7Y1.7_N256.*').name;
u17k = readmatrix(u_filename());
u_filename = dir('udata_X1.8Y1.8_N256.*').name;
u18k = readmatrix(u_filename());
u_filename = dir('udata_X1.9Y1.9_N256.*').name;
u19k = readmatrix(u_filename());
%un_filename = dir('u_norm_*').name;
%un15n = readmatrix(un15_filename());

tp = 150;
N = 512;
u11kv = NaN(N*N,tp);
un11k = NaN(1,tp);
u12kv = NaN(N*N,tp);
un12k = NaN(1,tp);
u13kv = NaN(N*N,tp);
un13k = NaN(1,tp);
u14kv = NaN(N*N,tp);
un14k = NaN(1,tp);
u15kv = NaN(N*N,tp);
un15k = NaN(1,tp);
u16kv = NaN(N*N,tp);
un16k = NaN(1,tp);
u17kv = NaN(N*N,tp);
un17k = NaN(1,tp);
u18kv = NaN(N*N,tp);
un18k = NaN(1,tp);
u19kv = NaN(N*N,tp);
un19k = NaN(1,tp);
for i = 1:tp
    u11kv(:,i) = reshape(u11k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un11k(1,i) = norm(u11kv(:,i))/N;
    u12kv(:,i) = reshape(u12k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un12k(1,i) = norm(u12kv(:,i))/N;
    u13kv(:,i) = reshape(u13k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un13k(1,i) = norm(u13kv(:,i))/N;
    u14kv(:,i) = reshape(u14k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un14k(1,i) = norm(u14kv(:,i))/N;
    u15kv(:,i) = reshape(u15k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un15k(1,i) = norm(u15kv(:,i))/N;
    u16kv(:,i) = reshape(u16k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un16k(1,i) = norm(u16kv(:,i))/N;
    u17kv(:,i) = reshape(u17k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un17k(1,i) = norm(u17kv(:,i))/N;
    u18kv(:,i) = reshape(u18k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un18k(1,i) = norm(u18kv(:,i))/N;
    u19kv(:,i) = reshape(u19k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un19k(1,i) = norm(u19kv(:,i))/N;
end

figure(1)
plot(dx(1, 1:101),normu11sinL(1:101, 1),'r','LineWidth',3)
hold on;
plot(dx(1, 1:101),normu12sinL(1:101, 1),'r','LineWidth',2)
plot(dx(1, 1:101),normu13sinL(1:101, 1),'r','LineWidth',1)
plot(dx(1, 1:101),normu14sinL(1:101, 1),'g','LineWidth',3)
plot(dx(1, 1:101),normu15sinL(1:101, 1),'g','LineWidth',2)
plot(dx(1, 1:101),normu16sinL(1:101, 1),'g','LineWidth',1)
plot(dx(1, 1:101),normu17sinL(1:101, 1),'b','LineWidth',3)
plot(dx(1, 1:101),normu18sinL(1:101, 1),'b','LineWidth',2)
plot(dx(1, 1:101),normu19sinL(1:101, 1),'b','LineWidth',1)
plot(dx(1, 1:101),normu11sinLmat(1:101, 1),'r','LineWidth',3)
hold on;
plot(dx(1, 1:101),normu12sinLmat(1:101, 1),'r','LineWidth',2)
plot(dx(1, 1:101),normu13sinLmat(1:101, 1),'r','LineWidth',1)
plot(dx(1, 1:101),normu14sinLmat(1:101, 1),'g','LineWidth',3)
plot(dx(1, 1:101),normu15sinLmat(1:101, 1),'g','LineWidth',2)
plot(dx(1, 1:101),normu16sinLmat(1:101, 1),'g','LineWidth',1)
plot(dx(1, 1:101),normu17sinLmat(1:101, 1),'b','LineWidth',3)
plot(dx(1, 1:101),normu18sinLmat(1:101, 1),'b','LineWidth',2)
plot(dx(1, 1:101),normu19sinLmat(1:101, 1),'b','LineWidth',1)
plot(un11k,'r','LineWidth',3)
hold on;
plot(un12k,'r','LineWidth',2)
plot(un13k,'r','LineWidth',1)
plot(un14k,'g','LineWidth',3)
plot(un15k,'g','LineWidth',2)
plot(un16k,'g','LineWidth',1)
plot(un17k,'b','LineWidth',3)
plot(un18k,'b','LineWidth',2)
plot(un19k,'b','LineWidth',1)
%}
set(gcf,'Position',[100 100 900 750])
        xlabel('$t$','Interpreter','latex'); 
        ylabel('$|| \phi(t) ||_{L^2}$','Interpreter','latex');

        fontsize(12,"points")
        set(gca,'fontsize', 14) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("2DKS growth rate for varying isotropic domains, IC sinL")
        legend("$\ell = 1.1$, F (N=512)", "$\ell = 1.2$, F (N=512)","$\ell = 1.3$, F (N=512)","$\ell = 1.4$, F (N=512)",...
            "$\ell = 1.5$, F (N=512)", "$\ell = 1.6$, F (N=512)","$\ell = 1.7$, F (N=512)","$\ell = 1.8$, F (N=512)",...
            "$\ell = 1.9$, F (N=512)", ...
            "$\ell = 1.1$, M (N=128)", "$\ell = 1.2$, M (N=128)","$\ell = 1.3$, M (N=128)","$\ell = 1.4$, M (N=128)",...
            "$\ell = 1.5$, M (N=128)", "$\ell = 1.6$, M (N=128)","$\ell = 1.7$, M (N=128)","$\ell = 1.8$, M (N=128)",...
            "$\ell = 1.9$, M (N=128)",...
            "$\ell = 1.1$, K (N=256)", "$\ell = 1.2$, K (N=256)","$\ell = 1.3$, K (N=256)","$\ell = 1.4$, K (N=256)",...
            "$\ell = 1.5$, K (N=256)", "$\ell = 1.6$, K (N=256)","$\ell = 1.7$, K (N=256)","$\ell = 1.8$, K (N=256)",...
            "$\ell = 1.9$, K (N=256)", ...
            'Location','southeast','NumColumns',3,'Interpreter','latex')

pts = 300;
diffL2_11kf = NaN(pts,1);
diffL2_12kf = NaN(pts,1);
diffL2_13kf = NaN(pts,1);
diffL2_14kf = NaN(pts,1);
diffL2_15kf = NaN(pts,1);
diffL2_16kf = NaN(pts,1);
diffL2_17kf = NaN(pts,1);
diffL2_18kf = NaN(pts,1);
diffL2_19kf = NaN(pts,1);
for j = 1:pts
    diffL2_11kf(j,1) = norm(normu11sinL(j+1,1)-un11k(1,j))/un11k(1,j);
    diffL2_12kf(j,1) = norm(normu12sinL(j+1,1)-un12k(1,j))/un12k(1,j);
    diffL2_13kf(j,1) = norm(normu13sinL(j+1,1)-un13k(1,j))/un13k(1,j);
    diffL2_14kf(j,1) = norm(normu14sinL(j+1,1)-un14k(1,j))/un14k(1,j);
    diffL2_15kf(j,1) = norm(normu15sinL(j+1,1)-un15k(1,j))/un15k(1,j);
    diffL2_16kf(j,1) = norm(normu16sinL(j+1,1)-un16k(1,j))/un16k(1,j);
    diffL2_17kf(j,1) = norm(normu17sinL(j+1,1)-un17k(1,j))/un17k(1,j);
    diffL2_18kf(j,1) = norm(normu18sinL(j+1,1)-un18k(1,j))/un18k(1,j);   
    diffL2_19kf(j,1) = norm(normu19sinL(j+1,1)-un19k(1,j))/un19k(1,j);
end
for j = 1
    if diffL2_11kf(j,1) == 0
        diffL2_11kf(j,1) = 1e-15;
    end
    if diffL2_12kf(j,1) == 0
        diffL2_12kf(j,1) = 1e-15;
    end
    if diffL2_13kf(j,1) == 0
        diffL2_13kf(j,1) = 1e-15;
    end
    if diffL2_14kf(j,1) == 0
        diffL2_14kf(j,1) = 1e-15;
    end
    if diffL2_15kf(j,1) == 0
        diffL2_15kf(j,1) = 1e-15;
    end
    if diffL2_16kf(j,1) == 0
        diffL2_16kf(j,1) = 1e-15;
    end
    if diffL2_17kf(j,1) == 0
        diffL2_17kf(j,1) = 1e-15;
    end    
    if diffL2_18kf(j,1) == 0
        diffL2_18kf(j,1) = 1e-15;
    end
    if diffL2_19kf(j,1) == 0
        diffL2_19kf(j,1) = 1e-15;
    end
end

figure(2)
semilogy(dx(1, 2:101),diffL2_11kf(1:150, 1),'r','LineWidth',2)
hold on;
semilogy(dx(1, 2:301),diffL2_12kf(1:300, 1),'r','LineWidth',2)
semilogy(dx(1, 2:101),diffL2_13kf(1:100, 1),'r','LineWidth',2)
semilogy(dx(1, 2:101),diffL2_14kf(1:100, 1),'g-o','LineWidth',2)
semilogy(dx(1, 2:101),diffL2_15kf(1:100, 1),'g-o','LineWidth',2)
semilogy(dx(1, 2:101),diffL2_16kf(1:100, 1),'g-o','LineWidth',2)
semilogy(dx(1, 2:101),diffL2_17kf(1:100, 1),'b-x','LineWidth',2)
semilogy(dx(1, 2:101),diffL2_18kf(1:100, 1),'b-x','LineWidth',2)
semilogy(dx(1, 2:101),diffL2_19kf(1:100, 1),'b-x','LineWidth',2)

set(gcf,'Position',[100 100 900 750])
        xlabel('$t$','Interpreter','latex'); 
        ylabel('$\frac{|| \phi_F(t) - \phi_K(t) ||_{L^2} }{ || \phi_K(t) ||_{L^2}}$','Interpreter','latex');

        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("Relative difference, KKP15 v. Fortran (N=512)")
        legend("$\ell = 1.1$", "$\ell = 1.2$", "$\ell = 1.3$", ...
            "$\ell = 1.4$", "$\ell = 1.5$", "$\ell = 1.6$", ...
            "$\ell = 1.7$", "$\ell = 1.8$", "$\ell = 1.9$", ...
                        'Location','northwest','NumColumns',3,'Interpreter','latex')

pts = 101;
diffL2_11mf = NaN(pts,1);
diffL2_12mf = NaN(pts,1);
diffL2_13mf = NaN(pts,1);
diffL2_14mf = NaN(pts,1);
diffL2_15mf = NaN(pts,1);
diffL2_16mf = NaN(pts,1);
diffL2_17mf = NaN(pts,1);
diffL2_18mf = NaN(pts,1);
diffL2_19mf = NaN(pts,1);
for j = 1:pts
    diffL2_11mf(j,1) = norm(normu11sinLmat(j,1)-normu11sinL(j,1))/normu11sinL(j,1);
    diffL2_12mf(j,1) = norm(normu12sinLmat(j,1)-normu12sinL(j,1))/normu12sinL(j,1);
    diffL2_13mf(j,1) = norm(normu13sinLmat(j,1)-normu13sinL(j,1))/normu13sinL(j,1);
    diffL2_14mf(j,1) = norm(normu14sinLmat(j,1)-normu14sinL(j,1))/normu14sinL(j,1);
    diffL2_15mf(j,1) = norm(normu15sinLmat(j,1)-normu15sinL(j,1))/normu15sinL(j,1);
    diffL2_16mf(j,1) = norm(normu16sinLmat(j,1)-normu16sinL(j,1))/normu16sinL(j,1);
    diffL2_17mf(j,1) = norm(normu17sinLmat(j,1)-normu17sinL(j,1))/normu17sinL(j,1);
    diffL2_18mf(j,1) = norm(normu18sinLmat(j,1)-normu18sinL(j,1))/normu18sinL(j,1);   
    diffL2_19mf(j,1) = norm(normu19sinLmat(j,1)-normu19sinL(j,1))/normu19sinL(j,1);
end
for j = 1
    if diffL2_11mf(j,1) == 0
        diffL2_11mf(j,1) = 1e-15;
    end
    if diffL2_12mf(j,1) == 0
        diffL2_12mf(j,1) = 1e-15;
    end
    if diffL2_13mf(j,1) == 0
        diffL2_13mf(j,1) = 1e-15;
    end
    if diffL2_14mf(j,1) == 0
        diffL2_14mf(j,1) = 1e-15;
    end
    if diffL2_15mf(j,1) == 0
        diffL2_15mf(j,1) = 1e-15;
    end
    if diffL2_16mf(j,1) == 0
        diffL2_16mf(j,1) = 1e-15;
    end
    if diffL2_17mf(j,1) == 0
        diffL2_17mf(j,1) = 1e-15;
    end    
    if diffL2_18mf(j,1) == 0
        diffL2_18mf(j,1) = 1e-15;
    end
    if diffL2_19mf(j,1) == 0
        diffL2_19mf(j,1) = 1e-15;
    end
end

figure(3)
semilogy(dx(1, 1:101),diffL2_11mf(1:101, 1),'r','LineWidth',2)
hold on;
semilogy(dx(1, 1:101),diffL2_12mf(1:101, 1),'r','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_13mf(1:101, 1),'r','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_14mf(1:101, 1),'g-o','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_15mf(1:101, 1),'g-o','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_16mf(1:101, 1),'g-o','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_17mf(1:101, 1),'b-x','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_18mf(1:101, 1),'b-x','LineWidth',2)
semilogy(dx(1, 1:101),diffL2_19mf(1:101, 1),'b-x','LineWidth',2)

set(gcf,'Position',[100 100 900 750])
        xlabel('$t$','Interpreter','latex'); 
        ylabel('$\frac{|| \phi_M(t) - \phi_F(t) ||_{L^2} }{ || \phi_F(t) ||_{L^2}}$','Interpreter','latex');

        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("Relative difference, Matlab (N=128) v. Fortran (N=512)")
        legend("$\ell = 1.1$", "$\ell = 1.2$", "$\ell = 1.3$", ...
            "$\ell = 1.4$", "$\ell = 1.5$", "$\ell = 1.6$", ...
            "$\ell = 1.7$", "$\ell = 1.8$", "$\ell = 1.9$", ...
                        'Location','northwest','NumColumns',3,'Interpreter','latex')

%% Test 5: sinL dynamics FK %%

%u140sinL = u_n;
%u150sinL = u_n;
%u160sinL = u_n;
%u170sinL = u_n;
%u180sinL = u_n;
%u190sinL = u_n;
%u198sinL = u_n;
%u210sinL = u_n;
%u225sinL = u_n;
%u235sinL = u_n;
%u250sinL = u_n;
%u267sinL = u_n;
%u280sinL = u_n;
%u301sinL = u_n;
%u320sinL = u_n;

% fortran

T = 300;
N = 512;
dt = 1e-3;

u_n = u140sinL;
ell1=1.4;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu140sinL = normu;

u_n = u150sinL;
ell1=1.5;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu150sinL = normu;

u_n = u160sinL;
ell1=1.6;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu160sinL = normu;

u_n = u170sinL;
ell1=1.7;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu170sinL = normu;

u_n = u180sinL;
ell1=1.8;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu180sinL = normu;

u_n = u190sinL;
ell1=1.9;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu190sinL = normu;

u_n = u198sinL;
ell1=1.98;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu198sinL = normu;

u_n = u210sinL;
ell1=2.1;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu210sinL = normu;

u_n = u225sinL;
ell1=2.25;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu225sinL = normu;

u_n = u235sinL;
ell1=2.35;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu235sinL = normu;

u_n = u250sinL;
ell1=2.50;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu250sinL = normu;

u_n = u267sinL;
ell1=2.67;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu267sinL = normu;

u_n = u280sinL;
ell1=2.80;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu280sinL = normu;

u_n = u301sinL;
ell1=3.01;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu301sinL = normu;

u_n = u320sinL;
ell1=3.20;
ell2 = ell1;
pts = size(u_n,2);
normu = NaN(pts,1);
for i = 1:pts
  normu(i,1) = norm(u_n(:,i))/N;
end
lognormu = log(normu);
dx = linspace(0,T,pts);
normu320sinL = normu;


u_filename = dir('udata_XY1.40.*').name;
u140k = readmatrix(u_filename());
u_filename = dir('udata_XY1.50.*').name;
u150k = readmatrix(u_filename());
u_filename = dir('udata_XY1.60.*').name;
u160k = readmatrix(u_filename());
u_filename = dir('udata_XY1.70.*').name;
u170k = readmatrix(u_filename());
u_filename = dir('udata_XY1.80.*').name;
u180k = readmatrix(u_filename());
u_filename = dir('udata_XY1.90.*').name;
u190k = readmatrix(u_filename());
u_filename = dir('udata_XY1.98.*').name;
u198k = readmatrix(u_filename());
u_filename = dir('udata_XY2.10.*').name;
u210k = readmatrix(u_filename());
u_filename = dir('udata_XY2.25.*').name;
u225k = readmatrix(u_filename());

u_filename = dir('udata_XY2.35.*').name;
u235k = readmatrix(u_filename());
u_filename = dir('udata_XY2.50.*').name;
u250k = readmatrix(u_filename());
u_filename = dir('udata_XY2.67.*').name;
u267k = readmatrix(u_filename());
u_filename = dir('udata_XY2.80.*').name;
u280k = readmatrix(u_filename());
u_filename = dir('udata_XY3.01.*').name;
u301k = readmatrix(u_filename());
u_filename = dir('udata_XY3.20.*').name;
u320k = readmatrix(u_filename());
%un_filename = dir('u_norm_*').name;
%un15n = readmatrix(un15_filename());

tp = 300;
N = 512;
u140kv = NaN(N*N,tp);
un140k = NaN(1,tp);
u150kv = NaN(N*N,tp);
un150k = NaN(1,tp);
u160kv = NaN(N*N,tp);
un160k = NaN(1,tp);
u170kv = NaN(N*N,tp);
un170k = NaN(1,tp);
u180kv = NaN(N*N,tp);
un180k = NaN(1,tp);
u190kv = NaN(N*N,tp);
un190k = NaN(1,tp);
u198kv = NaN(N*N,tp);
un198k = NaN(1,tp);
u210kv = NaN(N*N,tp);
un210k = NaN(1,tp);
u225kv = NaN(N*N,tp);
un225k = NaN(1,tp);

u235kv = NaN(N*N,tp);
un235k = NaN(1,tp);
u250kv = NaN(N*N,tp);
un250k = NaN(1,tp);
u267kv = NaN(N*N,tp);
un267k = NaN(1,tp);
u280kv = NaN(N*N,tp);
un280k = NaN(1,tp);
u301kv = NaN(N*N,tp);
un301k = NaN(1,tp);
u320kv = NaN(N*N,tp);
un320k = NaN(1,tp);
for i = 1:tp
    u140kv(:,i) = reshape(u140k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un140k(1,i) = norm(u140kv(:,i))/N;
    u150kv(:,i) = reshape(u150k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un150k(1,i) = norm(u150kv(:,i))/N;
    u160kv(:,i) = reshape(u160k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un160k(1,i) = norm(u160kv(:,i))/N;
    u170kv(:,i) = reshape(u170k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un170k(1,i) = norm(u170kv(:,i))/N;
    u180kv(:,i) = reshape(u180k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un180k(1,i) = norm(u180kv(:,i))/N;
    u190kv(:,i) = reshape(u190k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un190k(1,i) = norm(u190kv(:,i))/N;
    u198kv(:,i) = reshape(u198k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un198k(1,i) = norm(u198kv(:,i))/N;
    u210kv(:,i) = reshape(u210k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un210k(1,i) = norm(u210kv(:,i))/N;
    u225kv(:,i) = reshape(u225k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un225k(1,i) = norm(u225kv(:,i))/N;

    u235kv(:,i) = reshape(u235k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un235k(1,i) = norm(u235kv(:,i))/N;
    u250kv(:,i) = reshape(u250k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un250k(1,i) = norm(u250kv(:,i))/N;
    u267kv(:,i) = reshape(u267k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un267k(1,i) = norm(u267kv(:,i))/N;
    u280kv(:,i) = reshape(u280k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un280k(1,i) = norm(u280kv(:,i))/N;
    u301kv(:,i) = reshape(u301k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un301k(1,i) = norm(u301kv(:,i))/N;
    u320kv(:,i) = reshape(u320k((N*i)-(N-1):N*i,:) , [ N*N, 1 ] );
    un320k(1,i) = norm(u320kv(:,i))/N;
end

figure(1)
plot(dx(1, 1:251),normu140sinL(1:251, 1),'r','LineWidth',3)
hold on;
plot(dx(1, 1:251),normu150sinL(1:251, 1),'r','LineWidth',2)
plot(dx(1, 1:251),normu160sinL(1:251, 1),'r','LineWidth',1)
plot(dx(1, 1:251),normu170sinL(1:251, 1),'g','LineWidth',3)
plot(dx(1, 1:251),normu180sinL(1:251, 1),'g','LineWidth',2)
plot(dx(1, 1:251),normu190sinL(1:251, 1),'g','LineWidth',1)
plot(dx(1, 1:251),normu198sinL(1:251, 1),'b','LineWidth',3)
plot(dx(1, 1:251),normu210sinL(1:251, 1),'b','LineWidth',2)
plot(dx(1, 1:251),normu225sinL(1:251, 1),'b','LineWidth',1)

plot(dx(1, 1:251),normu235sinL(1:251, 1),'k','LineWidth',2)
plot(dx(1, 1:251),normu250sinL(1:251, 1),'m','LineWidth',2)
plot(dx(1, 1:251),normu267sinL(1:251, 1),'c','LineWidth',2)
plot(dx(1, 1:251),normu280sinL(1:251, 1),'k','LineWidth',1)
plot(dx(1, 1:251),normu301sinL(1:251, 1),'m','LineWidth',1)
plot(dx(1, 1:251),normu320sinL(1:251, 1),'c','LineWidth',1)

plot(un140k(1, 1:250),'r','LineWidth',3)
hold on;
plot(un150k(1, 1:250),'r','LineWidth',2)
plot(un160k(1, 1:250),'r','LineWidth',1)
plot(un170k(1, 1:250),'g','LineWidth',3)
plot(un180k(1, 1:250),'g','LineWidth',2)
plot(un190k(1, 1:250),'g','LineWidth',1)
plot(un198k(1, 1:250),'b','LineWidth',3)
plot(un210k(1, 1:250),'b','LineWidth',2)
plot(un225k(1, 1:250),'b','LineWidth',1)

plot(un235k(1, 1:250),'k','LineWidth',2)
plot(un250k(1, 1:250),'m','LineWidth',2)
plot(un267k(1, 1:250),'c','LineWidth',2)
plot(un280k(1, 1:250),'k','LineWidth',1)
plot(un301k(1, 1:250),'m','LineWidth',1)
plot(un320k(1, 1:250),'c','LineWidth',1)
%}
set(gcf,'Position',[100 100 900 750])
        xlabel('$t$','Interpreter','latex'); 
        ylabel('$|| \phi(t) ||_{L^2}$','Interpreter','latex');

        fontsize(12,"points")
        set(gca,'fontsize', 14) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("2DKS growth rate for varying isotropic domains, IC sinL")
        legend("$\ell = 1.40$, F", "$\ell = 1.50$, F","$\ell = 1.60$, F","$\ell = 1.70$, F",...
            "$\ell = 1.80$, F", "$\ell = 1.90$, F","$\ell = 1.98$, F","$\ell = 2.10$, F",...
            "$\ell = 2.25$, F", "$\ell = 2.35$, F","$\ell = 2.50$, F","$\ell = 2.67$, F",...
            "$\ell = 2.80$, F","$\ell = 3.01$, F","$\ell = 3.20$, F",...
            "$\ell = 1.40$, K", "$\ell = 1.50$, K","$\ell = 1.60$, K","$\ell = 1.70$, K",...
            "$\ell = 1.80$, K", "$\ell = 1.90$, K","$\ell = 1.98$, K","$\ell = 2.10$, K",...
            "$\ell = 2.25$, K", "$\ell = 2.35$, K","$\ell = 2.50$, K","$\ell = 2.67$, K",...
            "$\ell = 2.80$, K","$\ell = 3.01$, K","$\ell = 3.20$, K",...
            'Location','southoutside','NumColumns',5,'Interpreter','latex')

pts = 300;
diffL2_140kf = NaN(pts,1);
diffL2_150kf = NaN(pts,1);
diffL2_160kf = NaN(pts,1);
diffL2_170kf = NaN(pts,1);
diffL2_180kf = NaN(pts,1);
diffL2_190kf = NaN(pts,1);
diffL2_198kf = NaN(pts,1);
diffL2_210kf = NaN(pts,1);
diffL2_225kf = NaN(pts,1);

diffL2_235kf = NaN(pts,1);
diffL2_250kf = NaN(pts,1);
diffL2_267kf = NaN(pts,1);
diffL2_280kf = NaN(pts,1);
diffL2_301kf = NaN(pts,1);
diffL2_320kf = NaN(pts,1);
for j = 1:pts
    diffL2_140kf(j,1) = norm(normu140sinL(j+1,1)-un140k(1,j))/un140k(1,j);
    diffL2_150kf(j,1) = norm(normu150sinL(j+1,1)-un150k(1,j))/un150k(1,j);
    diffL2_160kf(j,1) = norm(normu160sinL(j+1,1)-un160k(1,j))/un160k(1,j);
    diffL2_170kf(j,1) = norm(normu170sinL(j+1,1)-un170k(1,j))/un170k(1,j);
    diffL2_180kf(j,1) = norm(normu180sinL(j+1,1)-un180k(1,j))/un180k(1,j);
    diffL2_190kf(j,1) = norm(normu190sinL(j+1,1)-un190k(1,j))/un190k(1,j);
    diffL2_198kf(j,1) = norm(normu198sinL(j+1,1)-un198k(1,j))/un198k(1,j);
    diffL2_210kf(j,1) = norm(normu210sinL(j+1,1)-un210k(1,j))/un210k(1,j);   
    diffL2_225kf(j,1) = norm(normu225sinL(j+1,1)-un225k(1,j))/un225k(1,j);

    diffL2_235kf(j,1) = norm(normu235sinL(j+1,1)-un235k(1,j))/un235k(1,j);
    diffL2_250kf(j,1) = norm(normu250sinL(j+1,1)-un250k(1,j))/un250k(1,j);
    diffL2_267kf(j,1) = norm(normu267sinL(j+1,1)-un267k(1,j))/un267k(1,j);
    diffL2_280kf(j,1) = norm(normu280sinL(j+1,1)-un280k(1,j))/un280k(1,j);
    diffL2_301kf(j,1) = norm(normu301sinL(j+1,1)-un301k(1,j))/un301k(1,j);
    diffL2_320kf(j,1) = norm(normu320sinL(j+1,1)-un320k(1,j))/un320k(1,j);
end
%{
for j = 1
    if diffL2_11kf(j,1) == 0
        diffL2_11kf(j,1) = 1e-15;
    end
    if diffL2_12kf(j,1) == 0
        diffL2_12kf(j,1) = 1e-15;
    end
    if diffL2_13kf(j,1) == 0
        diffL2_13kf(j,1) = 1e-15;
    end
    if diffL2_14kf(j,1) == 0
        diffL2_14kf(j,1) = 1e-15;
    end
    if diffL2_15kf(j,1) == 0
        diffL2_15kf(j,1) = 1e-15;
    end
    if diffL2_16kf(j,1) == 0
        diffL2_16kf(j,1) = 1e-15;
    end
    if diffL2_17kf(j,1) == 0
        diffL2_17kf(j,1) = 1e-15;
    end    
    if diffL2_18kf(j,1) == 0
        diffL2_18kf(j,1) = 1e-15;
    end
    if diffL2_19kf(j,1) == 0
        diffL2_19kf(j,1) = 1e-15;
    end
end
%}
figure(2)
semilogy(dx(1, 2:251),diffL2_140kf(1:250, 1),'r','LineWidth',2)
hold on;
semilogy(dx(1, 2:251),diffL2_150kf(1:250, 1),'g','LineWidth',2)
semilogy(dx(1, 2:251),diffL2_160kf(1:250, 1),'b','LineWidth',2)
semilogy(dx(1, 2:251),diffL2_170kf(1:250, 1),'r-o','LineWidth',2)
semilogy(dx(1, 2:251),diffL2_180kf(1:250, 1),'g-o','LineWidth',2)
semilogy(dx(1, 2:251),diffL2_190kf(1:250, 1),'b-o','LineWidth',2)
semilogy(dx(1, 2:251),diffL2_198kf(1:250, 1),'r-x','LineWidth',2)
semilogy(dx(1, 2:251),diffL2_210kf(1:250, 1),'g-x','LineWidth',2)
semilogy(dx(1, 2:251),diffL2_225kf(1:250, 1),'b-x','LineWidth',2)

semilogy(dx(1, 2:251),diffL2_235kf(1:250, 1),'k','LineWidth',2)
semilogy(dx(1, 2:251),diffL2_250kf(1:250, 1),'m','LineWidth',2)
semilogy(dx(1, 2:251),diffL2_267kf(1:250, 1),'c','LineWidth',2)
semilogy(dx(1, 2:251),diffL2_280kf(1:250, 1),'k-o','LineWidth',2)
semilogy(dx(1, 2:251),diffL2_301kf(1:250, 1),'m-o','LineWidth',2)
semilogy(dx(1, 2:251),diffL2_320kf(1:250, 1),'c-o','LineWidth',2)

set(gcf,'Position',[100 100 900 750])
        xlabel('$t$','Interpreter','latex'); 
        ylabel('$\frac{|| \phi_F(t) - \phi_K(t) ||_{L^2} }{ || \phi_K(t) ||_{L^2}}$','Interpreter','latex');

        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("Relative difference, KKP15 v. Fortran (N=512)")
        legend("$\ell = 1.40$", "$\ell = 1.50$","$\ell = 1.60$","$\ell = 1.70$",...
            "$\ell = 1.80$", "$\ell = 1.90$","$\ell = 1.98$","$\ell = 2.10$",...
            "$\ell = 2.25$", "$\ell = 2.35$","$\ell = 2.50$","$\ell = 2.67$",...
            "$\ell = 2.80$","$\ell = 3.01$","$\ell = 3.20$",...
                        'Location','southoutside','NumColumns',5,'Interpreter','latex')


%% GIF for KKP15 %%

dt = 1e-3;
T = 250;

ell1 = 1.4;
ell2 = ell1;
L1 = 2*pi*ell1;
L2 = 2*pi*ell2;
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

%u_n = u140kv;
u_n = u140sinL;

h = figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size

evolutiontiles_file = [pwd '/evolutiontiles_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(ell2) '_lY' num2str(ell1) '.gif'];
for i = 1 : tp + 1

    currentT = (i-1)/(tp)*T;
    subplot(1,2,1);
    % Draw surface plot
    %u_i = u_n((N*i)-(N-1):N*i,:);
    u_i = reshape( u_n(:,i) , [ N , N ] );
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );


    view([37.5,30]);
    drawnow

    subplot(1,2,2);
    % Draw surface plot
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );

    view(2);
    drawnow

    sgtitle(['2DKS DNS, x_1 = 2\pi(' num2str(ell1) '), x_2 = 2\pi(' num2str(ell2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])

    if i == 1
        %gif(evolution3D_file)
        gif(evolutiontiles_file)
    else
        gif
    end

end
close all

ell1 = 1.5;
ell2 = ell1;
L1 = 2*pi*ell1;
L2 = 2*pi*ell2;
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

%u_n = u140kv;
u_n = u150sinL;

h = figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size

evolutiontiles_file = [pwd '/evolutiontiles_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(ell2) '_lY' num2str(ell1) '.gif'];
for i = 1 : tp + 1

    currentT = (i-1)/(tp)*T;
    subplot(1,2,1);
    % Draw surface plot
    %u_i = u_n((N*i)-(N-1):N*i,:);
    u_i = reshape( u_n(:,i) , [ N , N ] );
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );


    view([37.5,30]);
    drawnow

    subplot(1,2,2);
    % Draw surface plot
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );

    view(2);
    drawnow

    sgtitle(['2DKS DNS, x_1 = 2\pi(' num2str(ell1) '), x_2 = 2\pi(' num2str(ell2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])

    if i == 1
        %gif(evolution3D_file)
        gif(evolutiontiles_file)
    else
        gif
    end

end
close all

ell1 = 1.6;
ell2 = ell1;
L1 = 2*pi*ell1;
L2 = 2*pi*ell2;
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

%u_n = u140kv;
u_n = u160sinL;

h = figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size

evolutiontiles_file = [pwd '/evolutiontiles_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(ell2) '_lY' num2str(ell1) '.gif'];
for i = 1 : tp + 1

    currentT = (i-1)/(tp)*T;
    subplot(1,2,1);
    % Draw surface plot
    %u_i = u_n((N*i)-(N-1):N*i,:);
    u_i = reshape( u_n(:,i) , [ N , N ] );
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );


    view([37.5,30]);
    drawnow

    subplot(1,2,2);
    % Draw surface plot
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );

    view(2);
    drawnow

    sgtitle(['2DKS DNS, x_1 = 2\pi(' num2str(ell1) '), x_2 = 2\pi(' num2str(ell2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])

    if i == 1
        %gif(evolution3D_file)
        gif(evolutiontiles_file)
    else
        gif
    end

end
close all

ell1 = 1.7;
ell2 = ell1;
L1 = 2*pi*ell1;
L2 = 2*pi*ell2;
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

%u_n = u140kv;
u_n = u170sinL;

h = figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size

evolutiontiles_file = [pwd '/evolutiontiles_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(ell2) '_lY' num2str(ell1) '.gif'];
for i = 1 : tp + 1

    currentT = (i-1)/(tp)*T;
    subplot(1,2,1);
    % Draw surface plot
    %u_i = u_n((N*i)-(N-1):N*i,:);
    u_i = reshape( u_n(:,i) , [ N , N ] );
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );


    view([37.5,30]);
    drawnow

    subplot(1,2,2);
    % Draw surface plot
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );

    view(2);
    drawnow

    sgtitle(['2DKS DNS, x_1 = 2\pi(' num2str(ell1) '), x_2 = 2\pi(' num2str(ell2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])

    if i == 1
        %gif(evolution3D_file)
        gif(evolutiontiles_file)
    else
        gif
    end

end
close all

ell1 = 1.8;
ell2 = ell1;
L1 = 2*pi*ell1;
L2 = 2*pi*ell2;
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

%u_n = u140kv;
u_n = u180sinL;

h = figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size

evolutiontiles_file = [pwd '/evolutiontiles_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(ell2) '_lY' num2str(ell1) '.gif'];
for i = 1 : tp + 1

    currentT = (i-1)/(tp)*T;
    subplot(1,2,1);
    % Draw surface plot
    %u_i = u_n((N*i)-(N-1):N*i,:);
    u_i = reshape( u_n(:,i) , [ N , N ] );
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );


    view([37.5,30]);
    drawnow

    subplot(1,2,2);
    % Draw surface plot
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );

    view(2);
    drawnow

    sgtitle(['2DKS DNS, x_1 = 2\pi(' num2str(ell1) '), x_2 = 2\pi(' num2str(ell2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])

    if i == 1
        %gif(evolution3D_file)
        gif(evolutiontiles_file)
    else
        gif
    end

end
close all

ell1 = 1.9;
ell2 = ell1;
L1 = 2*pi*ell1;
L2 = 2*pi*ell2;
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

%u_n = u140kv;
u_n = u190sinL;

h = figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size

evolutiontiles_file = [pwd '/evolutiontiles_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(ell2) '_lY' num2str(ell1) '.gif'];
for i = 1 : tp + 1

    currentT = (i-1)/(tp)*T;
    subplot(1,2,1);
    % Draw surface plot
    %u_i = u_n((N*i)-(N-1):N*i,:);
    u_i = reshape( u_n(:,i) , [ N , N ] );
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );


    view([37.5,30]);
    drawnow

    subplot(1,2,2);
    % Draw surface plot
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );

    view(2);
    drawnow

    sgtitle(['2DKS DNS, x_1 = 2\pi(' num2str(ell1) '), x_2 = 2\pi(' num2str(ell2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])

    if i == 1
        %gif(evolution3D_file)
        gif(evolutiontiles_file)
    else
        gif
    end

end
close all

ell1 = 1.98;
ell2 = ell1;
L1 = 2*pi*ell1;
L2 = 2*pi*ell2;
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

%u_n = u140kv;
u_n = u198sinL;

h = figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size

evolutiontiles_file = [pwd '/evolutiontiles_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(ell2) '_lY' num2str(ell1) '.gif'];
for i = 1 : tp + 1

    currentT = (i-1)/(tp)*T;
    subplot(1,2,1);
    % Draw surface plot
    %u_i = u_n((N*i)-(N-1):N*i,:);
    u_i = reshape( u_n(:,i) , [ N , N ] );
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );


    view([37.5,30]);
    drawnow

    subplot(1,2,2);
    % Draw surface plot
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );

    view(2);
    drawnow

    sgtitle(['2DKS DNS, x_1 = 2\pi(' num2str(ell1) '), x_2 = 2\pi(' num2str(ell2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])

    if i == 1
        %gif(evolution3D_file)
        gif(evolutiontiles_file)
    else
        gif
    end

end
close all

ell1 = 2.1;
ell2 = ell1;
L1 = 2*pi*ell1;
L2 = 2*pi*ell2;
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

%u_n = u140kv;
u_n = u210sinL;

h = figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size

evolutiontiles_file = [pwd '/evolutiontiles_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(ell2) '_lY' num2str(ell1) '.gif'];
for i = 1 : tp + 1

    currentT = (i-1)/(tp)*T;
    subplot(1,2,1);
    % Draw surface plot
    %u_i = u_n((N*i)-(N-1):N*i,:);
    u_i = reshape( u_n(:,i) , [ N , N ] );
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );


    view([37.5,30]);
    drawnow

    subplot(1,2,2);
    % Draw surface plot
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );

    view(2);
    drawnow

    sgtitle(['2DKS DNS, x_1 = 2\pi(' num2str(ell1) '), x_2 = 2\pi(' num2str(ell2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])

    if i == 1
        %gif(evolution3D_file)
        gif(evolutiontiles_file)
    else
        gif
    end

end
close all

ell1 = 2.25;
ell2 = ell1;
L1 = 2*pi*ell1;
L2 = 2*pi*ell2;
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

%u_n = u140kv;
u_n = u225sinL;

h = figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size

evolutiontiles_file = [pwd '/evolutiontiles_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(ell2) '_lY' num2str(ell1) '.gif'];
for i = 1 : tp + 1

    currentT = (i-1)/(tp)*T;
    subplot(1,2,1);
    % Draw surface plot
    %u_i = u_n((N*i)-(N-1):N*i,:);
    u_i = reshape( u_n(:,i) , [ N , N ] );
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );


    view([37.5,30]);
    drawnow

    subplot(1,2,2);
    % Draw surface plot
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );

    view(2);
    drawnow

    sgtitle(['2DKS DNS, x_1 = 2\pi(' num2str(ell1) '), x_2 = 2\pi(' num2str(ell2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])

    if i == 1
        %gif(evolution3D_file)
        gif(evolutiontiles_file)
    else
        gif
    end

end
close all

ell1 = 2.35;
ell2 = ell1;
L1 = 2*pi*ell1;
L2 = 2*pi*ell2;
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

%u_n = u140kv;
u_n = u235sinL;

h = figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size

evolutiontiles_file = [pwd '/evolutiontiles_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(ell2) '_lY' num2str(ell1) '.gif'];
for i = 1 : tp + 1

    currentT = (i-1)/(tp)*T;
    subplot(1,2,1);
    % Draw surface plot
    %u_i = u_n((N*i)-(N-1):N*i,:);
    u_i = reshape( u_n(:,i) , [ N , N ] );
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );


    view([37.5,30]);
    drawnow

    subplot(1,2,2);
    % Draw surface plot
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );

    view(2);
    drawnow

    sgtitle(['2DKS DNS, x_1 = 2\pi(' num2str(ell1) '), x_2 = 2\pi(' num2str(ell2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])

    if i == 1
        %gif(evolution3D_file)
        gif(evolutiontiles_file)
    else
        gif
    end

end
close all

ell1 = 2.5;
ell2 = ell1;
L1 = 2*pi*ell1;
L2 = 2*pi*ell2;
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

%u_n = u140kv;
u_n = u250sinL;

h = figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size

evolutiontiles_file = [pwd '/evolutiontiles_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(ell2) '_lY' num2str(ell1) '.gif'];
for i = 1 : tp + 1

    currentT = (i-1)/(tp)*T;
    subplot(1,2,1);
    % Draw surface plot
    %u_i = u_n((N*i)-(N-1):N*i,:);
    u_i = reshape( u_n(:,i) , [ N , N ] );
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );


    view([37.5,30]);
    drawnow

    subplot(1,2,2);
    % Draw surface plot
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );

    view(2);
    drawnow

    sgtitle(['2DKS DNS, x_1 = 2\pi(' num2str(ell1) '), x_2 = 2\pi(' num2str(ell2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])

    if i == 1
        %gif(evolution3D_file)
        gif(evolutiontiles_file)
    else
        gif
    end

end
close all

ell1 = 2.67;
ell2 = ell1;
L1 = 2*pi*ell1;
L2 = 2*pi*ell2;
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

%u_n = u140kv;
u_n = u267sinL;

h = figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size

evolutiontiles_file = [pwd '/evolutiontiles_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(ell2) '_lY' num2str(ell1) '.gif'];
for i = 1 : tp + 1

    currentT = (i-1)/(tp)*T;
    subplot(1,2,1);
    % Draw surface plot
    %u_i = u_n((N*i)-(N-1):N*i,:);
    u_i = reshape( u_n(:,i) , [ N , N ] );
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );


    view([37.5,30]);
    drawnow

    subplot(1,2,2);
    % Draw surface plot
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );

    view(2);
    drawnow

    sgtitle(['2DKS DNS, x_1 = 2\pi(' num2str(ell1) '), x_2 = 2\pi(' num2str(ell2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])

    if i == 1
        %gif(evolution3D_file)
        gif(evolutiontiles_file)
    else
        gif
    end

end
close all

ell1 = 2.8;
ell2 = ell1;
L1 = 2*pi*ell1;
L2 = 2*pi*ell2;
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

%u_n = u140kv;
u_n = u280sinL;

h = figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size

evolutiontiles_file = [pwd '/evolutiontiles_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(ell2) '_lY' num2str(ell1) '.gif'];
for i = 1 : tp + 1

    currentT = (i-1)/(tp)*T;
    subplot(1,2,1);
    % Draw surface plot
    %u_i = u_n((N*i)-(N-1):N*i,:);
    u_i = reshape( u_n(:,i) , [ N , N ] );
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );


    view([37.5,30]);
    drawnow

    subplot(1,2,2);
    % Draw surface plot
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );

    view(2);
    drawnow

    sgtitle(['2DKS DNS, x_1 = 2\pi(' num2str(ell1) '), x_2 = 2\pi(' num2str(ell2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])

    if i == 1
        %gif(evolution3D_file)
        gif(evolutiontiles_file)
    else
        gif
    end

end
close all

ell1 = 3.01;
ell2 = ell1;
L1 = 2*pi*ell1;
L2 = 2*pi*ell2;
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

%u_n = u140kv;
u_n = u301sinL;

h = figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size

evolutiontiles_file = [pwd '/evolutiontiles_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(ell2) '_lY' num2str(ell1) '.gif'];
for i = 1 : tp + 1

    currentT = (i-1)/(tp)*T;
    subplot(1,2,1);
    % Draw surface plot
    %u_i = u_n((N*i)-(N-1):N*i,:);
    u_i = reshape( u_n(:,i) , [ N , N ] );
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );


    view([37.5,30]);
    drawnow

    subplot(1,2,2);
    % Draw surface plot
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );

    view(2);
    drawnow

    sgtitle(['2DKS DNS, x_1 = 2\pi(' num2str(ell1) '), x_2 = 2\pi(' num2str(ell2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])

    if i == 1
        %gif(evolution3D_file)
        gif(evolutiontiles_file)
    else
        gif
    end

end
close all

ell1 = 3.2;
ell2 = ell1;
L1 = 2*pi*ell1;
L2 = 2*pi*ell2;
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

%u_n = u140kv;
u_n = u320sinL;

h = figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size

evolutiontiles_file = [pwd '/evolutiontiles_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(ell2) '_lY' num2str(ell1) '.gif'];
for i = 1 : tp + 1

    currentT = (i-1)/(tp)*T;
    subplot(1,2,1);
    % Draw surface plot
    %u_i = u_n((N*i)-(N-1):N*i,:);
    u_i = reshape( u_n(:,i) , [ N , N ] );
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );


    view([37.5,30]);
    drawnow

    subplot(1,2,2);
    % Draw surface plot
    surfc(x1,x2,u_i);
    xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
    shading(gca,'interp')
    set(gcf,'color','white')
    set(gca,'color','white')
    colormap(redblue)
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );

    view(2);
    drawnow

    sgtitle(['2DKS DNS, x_1 = 2\pi(' num2str(ell1) '), x_2 = 2\pi(' num2str(ell2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])

    if i == 1
        %gif(evolution3D_file)
        gif(evolutiontiles_file)
    else
        gif
    end

end
close all

%% Theoretical growth rates %%%

lampts = 1e4 + 1;
Lmax = 6.0;
x_lam = linspace(1,Lmax,lampts);

% isotropic case
lam_max_iso = NaN(lampts,floor(Lmax) + 1);
lam_mode_iso = zeros(lampts,floor(Lmax+1)^2);
lam_max_i = NaN(lampts,1);
wavevector = NaN(floor(Lmax+1)^2,3);
for i = 1:lampts
    ell = x_lam(i);
    lamk = NaN(floor(ell)+1,floor(ell)+1);
    push = 0;
    for k1 = 0:floor(ell)
        index = 1;
        for k2 = 0:floor(ell)
            wavevector(push + index,1:2) = [ k1, k2 ] ;
            wavevector(push + index,3) = k1^2 + k2^2  ;
            lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
            lam_mode_iso(i, push + index ) = lamk(k1+1,k2+1);
            index = index + 1;
        end
        push = push + floor(Lmax+1);
    end
    lam_max_i(i,1) = max(max(lamk));
    for j = floor(ell) + 1
        lam_max_iso(i,1:j) = maxk(lamk(:),j);
    end
end

% decay
lampts = 1e4 + 1;
Lmax = 1.0;
x_lam = linspace(0,Lmax,lampts);
lam_max_d = NaN(lampts,1);
for i = 1:lampts
    ell = x_lam(i);
    lamk = NaN(5,5);
    for k1 = 0:5
        for k2 = 0:5
            lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
        end
    end
    lamk(1,1) = -inf;
    lam_max_d(i,1) = max(max(lamk));
end

lam_max_1 = NaN(lampts,1);
ell1 = 1; % anisotropic with 1 growing mode in dimension 1
for i = 1:lampts
    ell = x_lam(i);
    lamk = NaN(floor(ell),floor(ell));
    for k1 = 0:floor(ell)
        for k2 = 0:floor(ell)
            lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
        end
    end
    lam_max_1(i,1) = max(max(lamk));
end

lam_max_2 = NaN(lampts,1);
ell1 = 2; % anisotropic with 2 growing modes in dimension 1
for i = 1:lampts
    ell = x_lam(i);
    lamk = NaN(floor(ell),floor(ell));
    for k1 = 0:floor(ell)
        for k2 = 0:floor(ell)
            lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
        end
    end
    lam_max_2(i,1) = max(max(lamk));
end

lam_max_3 = NaN(lampts,1);
ell1 = 3; % anisotropic with 3 growing modes in dimension 1
for i = 1:lampts
    ell = x_lam(i);
    lamk = NaN(floor(ell),floor(ell));
    for k1 = 0:floor(ell)
        for k2 = 0:floor(ell)
            lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
        end
    end
    lam_max_3(i,1) = max(max(lamk));
end

lam_max_4 = NaN(lampts,1);
ell1 = 4; % anisotropic with 4 growing modes in dimension 1
for i = 1:lampts
    ell = x_lam(i);
    lamk = NaN(floor(ell),floor(ell));
    for k1 = 0:floor(ell)
        for k2 = 0:floor(ell)
            lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
        end
    end
    lam_max_4(i,1) = max(max(lamk));
end

lam_max_5 = NaN(lampts,1);
ell1 = 5; % anisotropic with 5 growing modes in dimension 1
for i = 1:lampts
    ell = x_lam(i);
    lamk = NaN(floor(ell),floor(ell));
    for k1 = 0:floor(ell)
        for k2 = 0:floor(ell)
            lamk(k1+1,k2+1) = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;
        end
    end
    lam_max_5(i,1) = max(max(lamk));
end

figure(3)
plot(x_lam,lam_max_i,'r','LineWidth',2)
hold on;
plot(x_lam,lam_max_2,'b','LineWidth',2)
plot(x_lam,lam_max_3,'m','LineWidth',2)
plot(x_lam,lam_max_4,'c','LineWidth',2)
plot(x_lam,lam_max_5,'k','LineWidth',2)
set(gcf,'Position',[100 100 900 750])
xlabel('Domain factor $\ell_2$','Interpreter','latex'); 
ylabel('Decay rate $\lambda^*$','Interpreter','latex');
xlim([1 Lmax])
ylim([0.22 0.25])

fontsize(12,"points")
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')    
title("Solution decay rate v. Domain size")
legend("$\ell_1 = \ell_2$", "$\ell_1 = 1$", "$\ell_1 = 2$", ...
    "$\ell_1 = 3$", "$\ell_1 = 4$", "$\ell_1 = 5$",  ...
                'Location','southeast','NumColumns',1,'Interpreter','latex')

figure(4)
plot(x_lam,lam_max_iso(:,1),'r','LineWidth',3)
hold on;
plot(x_lam,lam_max_iso(:,2),'b','LineWidth',3)
plot(x_lam,lam_max_iso(:,3),'m','LineWidth',3)
plot(x_lam,lam_max_iso(:,4),'c','LineWidth',2)
plot(x_lam,lam_max_iso(:,5),'k','LineWidth',2)
plot(x_lam,lam_max_iso(:,6),'g','LineWidth',2)
%plot(x_lam,lam_max_iso(:,7),'r','LineWidth',1)
%plot(x_lam,lam_max_iso(:,8),'b','LineWidth',1)
%plot(x_lam,lam_max_iso(:,9),'m','LineWidth',1)
set(gcf,'Position',[100 100 900 750])
xlabel('Domain factor $\ell_2$','Interpreter','latex'); 
xlim([1 Lmax])
ylim([0 0.25])
ylabel('Growth rate $\lambda^*$','Interpreter','latex');

fontsize(12,"points")
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')    
title("Solution growth rate v. Isotropic domain size")
legend("$\lambda_{1}$", "$\lambda_{2}$", ...
    "$\lambda_{3}$", "$\lambda_{4}$", "$\lambda_{5}$", "$\lambda_{6}$", ...
    "$\lambda_{7}$", "$\lambda_{8}$", "$\lambda_{9}$",  ...
                'Location','southoutside','NumColumns',9,'Interpreter','latex')

figure(5)
plot(x_lam,lam_mode_iso(:,1),'LineWidth',0.5)
deja = NaN(size(wavevector,1),1);
hold on;
for i = 2:length(lamk)^2
    plot(x_lam,lam_mode_iso(:,i),'.','LineWidth',0.5)
    [ val, ind ] = max(lam_mode_iso(:,i));
    wavec = ['(' num2str(wavevector(i,1)) ',' num2str(wavevector(i,2)) ')']; 
    if val > 0.24 && x_lam(ind) < 5.7 && ~ismember(wavevector(i,3), deja)
        %text(x_lam(round(ind))-0.125,lam_mode_iso(round(ind),i)+.005, wavec );
        text(x_lam(round(ind))-0.125,lam_mode_iso(round(ind),i)+.00075, wavec );
        deja(i) = wavevector(i,3);
    end
end
plot(x_lam,lam_max_i,'rx','LineWidth',1)
set(gcf,'Position',[100 100 900 750])
xlabel('Domain factor $\ell$','Interpreter','latex'); 
xlim([1.0 Lmax])
%ylim([0.001 0.26])
ylim([0.22 0.2515])
ylabel('Growth rate $\lambda$','Interpreter','latex');

%xline(sqrt(2),'b--','$\ell=\ell_{\max}$', 'Interpreter','latex','LabelOrientation','horizontal','LabelVerticalAlignment','bottom')
%xline(sqrt(3),'k--','$\ell=\ell_c$', 'Interpreter','latex','LabelOrientation','horizontal','LabelVerticalAlignment','bottom')
xline(sqrt(2),'b--','LineWidth',1)
xline(sqrt(3),'k--','LineWidth',1)

fontsize(12,"points")
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')    
title("Growing modes in varying isotropic domain sizes",'Interpreter','latex')
%legend("Fourier mode", 'Location','southeast','NumColumns',9,'Interpreter','latex')
