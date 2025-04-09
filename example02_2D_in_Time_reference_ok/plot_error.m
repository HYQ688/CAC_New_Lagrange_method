close all;
clear; clc;

domain.left   = 0;
domain.right  = 2*pi;
domain.bottom = 0;
domain.top    = 2*pi;

Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;

N  = 64;
Nx = N;
Ny = N;

Plot = 1;

if Plot == 1 

    % scheme1 = 'linear';
    scheme1 = 'nonlinear';

%     scheme2 = '_1st';   % First-order scheme
%     scheme2 = '_2cn';   % Second-order CN scheme
    scheme2 = '_bdf2';  % Second-order BDF scheme

scheme = [scheme1 , scheme2];

else
    % sch1 = 'linear';
    sch1 = 'nonlinear';

    % sch2 = '_1st'; 
    sch2 = '_bdf2';

%     sch3 = 'linear';
    sch3 = 'nonlinear';

    % sch4 = '_1st'; 
    sch4 = '_bdf2';

    scheme1 = [sch1 , sch2];
    scheme2 = [sch3 , sch4];
end


% Parameters
para.epsilon = 6*pi/128;
para.gamma = 1;
para.Re = 1;
para.lambda = 0.01;
para.S1 = 4;
para.S2 = 4;
para.S3 = 1;
PDE = 'data1';

if 1 == strcmp(PDE,'data1')
    T = 2;
    dt_array = 0.01./2.^(1:6)';
    dt_array1 = 0.003./2.^(1:6)';
    dt_ref = 1e-6;
    para.C0 = 100;
    pde = ex03_1_Vesicles_data(para);
end

if 1 == strcmp(PDE,'data2')
    T = 0.2;
    dt_array = 0.01./2.^(0:6)';
    dt_ref = 1e-6;
    para.C0 = 100;
    pde = ex03_2_Vesicles_data(para);
end

% Time: dt T

t0 = 0;
tsave = 0.2*T;

maxIt = length(dt_array);

%% Compute order of convergence(1)
if Plot == 1
error=zeros(maxIt,1);
order=zeros(maxIt,1);
if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
    name=[scheme,'_phi_e',num2str(pde.epsilon),'gamma=',num2str(pde.gamma),'S1=',num2str(pde.S1),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt_ref)];
    filename=[name '.mat'];
    load(filename,'phi','u','v','p');
    phi_exact = phi;
    u_exact = u;
    v_exact = v;
    p_exact = p;
    clear phi u v p;
else
    Lx = domain.right - domain.left;
    Ly = domain.top   - domain.bottom;
    hx = Lx/Nx;
    hy = Ly/Ny;
    x  = domain.left   + hx*(0:Nx-1);
    y  = domain.bottom + hy*(0:Ny-1);
    [xx,yy] = ndgrid(x,y);
    phi_exact = pde.exactphi(xx,yy,T);
    u_exact = pde.exactu(xx,yy,T);
    v_exact = pde.exactv(xx,yy,T);
    p_exact = pde.exactp(xx,yy,T);
end
for k = 1:maxIt
    dt = dt_array(k);
    name=[scheme,'_phi_e',num2str(pde.epsilon),'gamma=',num2str(pde.gamma),'S1=',num2str(pde.S1),...
        'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];
    filenamek=[name '.mat'];
    load(filenamek,'phi','u','v','p','hx','hy');
    err_phi    = fft2((phi_exact - phi).^2);
    error_phi(k,1) = sqrt(err_phi(1,1)*hx*hy);   % L2

    err_u = fft2((u_exact - u).^2 + (v_exact - v).^2);
    error_u(k,1) = sqrt(err_u(1,1)*hx*hy);   % L2

    err_p = fft2((p_exact - p).^2);
    error_p(k,1) = sqrt(err_p(1,1)*hx*hy);   % L2

    clear phi u v p;
end
order_phi(2:maxIt) = log(error_phi(1:maxIt-1)./error_phi(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));
order_u(2:maxIt)   = log(error_u(1:maxIt-1)./error_u(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));
order_p(2:maxIt)   = log(error_p(1:maxIt-1)./error_p(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));

%% Plot(1)
% figure(2)
hh=loglog(dt_array,error_phi,'*-','LineWidth',3,'MarkerSize',15);
hold on;
hh=loglog(dt_array,error_u, 'p-', 'LineWidth',3,'MarkerSize',15);
hh=loglog(dt_array,error_p, 'b-o','LineWidth',3,'MarkerSize',15);
hold on 
tf = strcmp(scheme2,'_1st');
if tf == 1
    hh=loglog(dt_array,dt_array1,'k-.','LineWidth',3,'MarkerSize',10);
    hold on
    h = legend('SAV LM scheme: $\phi$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
else 
    hh=loglog(dt_array,dt_array1.^2,'k-.','LineWidth',3,'MarkerSize',15);
    hold on
    h = legend('$Scheme2:\phi$','$Scheme2:\mathbf{u}$','$Scheme2:p$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
    % h = legend('$Scheme1:\phi$','$Scheme1:\mathbf{u}$','$Scheme1:p$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
end
set(h,'interpreter','latex','FontSize',15);
xlabel('Time step $\delta t$','Interpreter','latex','FontSize',18);
%  ylabel('Error\_max')
ylabel('$L^2$ error','Interpreter','latex','FontSize',18);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
% ylim([1e-8 5e-1])
ylim([1e-12 5e-3])
grid on;
hold on;

figname1 = ['D:\paper\phase-field\CAC_Vesicles_surface_LagrangeMultiplier_SAV\fig\',scheme,'_accuracy','.png'];
print(figname1,'-dpng', '-r300')

%% Save error and order(1)
name=[scheme,'_phi_e',num2str(pde.epsilon),'gamma=',num2str(pde.gamma),'S1=',num2str(pde.S1),...
        'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];
fileID = fopen([name,'.txt'],'w');
% fprintf(fileID,'%6s\n','%% Results');
% fprintf(fileID,'%6s\n','% dt	   &   Error_L2	   &  Order');
% A = [dt_array error];
% fprintf(fileID,'%.12f   %.4e   \n',A');
fprintf(fileID,'%.12f     %.4e      %.2f \n',[dt_array,error,order]');
fclose(fileID);

else 
%% Compute order of convergence(2)
error1=zeros(maxIt,1);
order1=zeros(maxIt,1);
error2=zeros(maxIt,1);
order2=zeros(maxIt,1);
if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
    name1=[scheme1,'_phi_e',num2str(pde.epsilon),'gamma',num2str(pde.gamma),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt_ref)];
    name2=[scheme2,'_phi_e',num2str(pde.epsilon),'gamma',num2str(pde.gamma),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt_ref)];
    filename1=[name1 '.mat'];
    filename2=[name2 '.mat'];
    load(filename1,'phi');
    phi_exact1 = phi;
    clear phi
    load(filename2,'phi');
    phi_exact2 = phi;
    clear phi;
else
    Lx = domain.right - domain.left;
    Ly = domain.top   - domain.bottom;
    hx = Lx/Nx;
    hy = Ly/Ny;
    x  = domain.left   + hx*(0:Nx-1);
    y  = domain.bottom + hy*(0:Ny-1);
    [xx,yy] = ndgrid(x,y);
    phi_exact = pde.exact(xx,yy,T);
end
for k = 1:maxIt
    dt = dt_array(k);
    name1=[scheme1,'_phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
        'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];
    name2=[scheme2,'_phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
        'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];
    filenamek1=[name1 '.mat'];
    load(filenamek1,'phi','hx','hy');
    err1    = fft2((phi_exact1 - phi).^2);
    error1(k,1) = sqrt(err1(1,1)*hx*hy);   % L2
    clear phi;
    filenamek2=[name2 '.mat'];
    load(filenamek2,'phi','hx','hy');
    err2    = fft2((phi_exact2 - phi).^2);
    error2(k,1) = sqrt(err2(1,1)*hx*hy);   % L2
    clear phi;
end
order1(2:maxIt) = log(error1(1:maxIt-1)./error1(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));
order2(2:maxIt) = log(error2(1:maxIt-1)./error2(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));

%% Plot(2)
figure(2)
hh=loglog(dt_array,error1,'*-','LineWidth',2,'MarkerSize',10);
hold on
hh=loglog(dt_array,error2,'p-','LineWidth',2,'MarkerSize',10);
hold on 
tf1 = strcmp(sch2,'_1st');
tf2 = strcmp(sch4,'_1st');
if tf1 == 1 && tf2 == 1
    hh=loglog(dt_array,dt_array,'k-.','LineWidth',3,'MarkerSize',10);
    hold on
    h = legend('First scheme: $\phi$','Second scheme: $\phi$','$\mathcal{O}(\delta t)$','Location','SouthEast');
elseif tf1 == 0 && tf2 == 0
    hh=loglog(dt_array,dt_array.^2,'k-.','LineWidth',3,'MarkerSize',10);
    hold on 
    h = legend('First scheme: $\phi$','Second scheme: $\phi$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
else 
    hh=loglog(dt_array,dt_array,'k-.','LineWidth',3,'MarkerSize',10);
    hold on
    hh=loglog(dt_array,dt_array.^2,'k-.','LineWidth',3,'MarkerSize',10);
    hold on
    h = legend('First scheme: $\phi$','Second scheme: $\phi$','$\mathcal{O}(\delta t)$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
end
% hh=loglog(dt_array,dt_array.^2,'k-.','LineWidth',3,'MarkerSize',10);
set(h,'interpreter','latex','FontSize',24);
xlabel('Time step $\delta t$','Interpreter','latex');
%  ylabel('Error\_max')
ylabel('$L^2$ error','Interpreter','latex');
xlim([2e-7 2e-4])
% ylim([1e-11 5e-0])
% ylim([1e-9 5e-5])
ylim([0 5e-3])
grid on;
hold on;

set(gca,'FontSize',24);
set(gca,'linewidth',1.5)

scheme = [scheme1, scheme2];

figname1 = ['D:\paper\phase-field\CAC_Vesicles_surface_LagrangeMultiplier_SAV\fig\',scheme1,'and',scheme2,'_accuracy','.png'];
% print(figname1,'-dpng', '-r300')

end