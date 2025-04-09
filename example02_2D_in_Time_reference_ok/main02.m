close all;
clear; clc;

% add path
addpath('../','-begin');

% Space: Domain and N
domain.left   = 0;
domain.right  = 2*pi;
domain.bottom = 0;
domain.top    = 2*pi;

Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;

N  = 128 ;
Nx = N;
Ny = N;

scheme1 = 'linear';
% scheme1 = 'nonlinear';

scheme2 = '_1st';   % First-order scheme
% scheme2 = '_2cn';   % Second-order CN scheme
% scheme2 = '_bdf2';  % Second-order BDF scheme

scheme = [scheme1 , scheme2];

% Parameters
para.epsilon = 6*pi/128;
% para.epsilon = 0.08;
para.gamma = 1;
para.Re = 1;
% para.lambda = 0.001;
para.lambda = 0.01;
para.S1 = 4;
para.S2 = 4;
para.S3 = 1;

PDE = 'data1';

if 1 == strcmp(PDE,'data1')
    T = 2;
    dt_array = 0.01./2.^(1:6)';
    dt_ref = 1e-6;
    para.C0 = 100;
    pde = ex03_1_Vesicles_data(para);
end

if 1 == strcmp(PDE,'data2')
    T = 0.2;
    dt_array = 0.01./2.^(0:7)';
    dt_ref = 1e-5;
    para.C0 = 100;
    pde = ex03_2_Vesicles_data(para);
end

% Time: dt T

t0 = 0;
tsave = 0.2*T;

maxIt = length(dt_array);

%% option
option.scheme = scheme;
option.plotflag  = 0;
option.printflag = 1;
option.vtkflag  = 0;
option.saveflag  = 0;
option.savefinal  = 1;
option.energyflag = 0;
option.tol = 1e-14;
option.tolit = 1e-11;
option.maxit = 2000;

if 1 == strcmp(scheme,'linear_1st')
    solver_fun = @CAC_Vesicle_with_NS_2D_linear_newLM_p_SAV_1st;
elseif 1 == strcmp(scheme,'linear_2cn')
    solver_fun = @CAC_Vesicle_2D_linear_newLM_p_SAV_2cn;
elseif 1 == strcmp(scheme,'linear_bdf2')
    solver_fun = @CAC_Vesicle_with_NS_2D_linear_newLM_p_SAV_bdf2;
elseif 1 == strcmp(scheme,'nonlinear_1st')
    solver_fun = @CAC_Vesicle_with_NS_2D_newLM_p_SAV_1st;
elseif 1 == strcmp(scheme,'nonlinear_2cn')
    solver_fun = @CAC_Vesicle_2D_newLM_p_SAV_2cn;
elseif 1 == strcmp(scheme,'nonlinear_bdf2')
    solver_fun = @CAC_Vesicle_with_NS_2D_newLM_p_SAV_bdf2;
end


%% Run:
% delete *.mat
if ~isfield(pde,'exactphi') || ~isfield(pde,'rhs1')
    time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
    solver_fun(pde,domain,Nx,Ny,time,option);
end
for k = 1:maxIt
    dt = dt_array(k);
    time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
    v2 = solver_fun(pde,domain,Nx,Ny,time,option);
end



%% Compute order of convergence
error_phi=zeros(maxIt,1);
order_phi=zeros(maxIt,1);
error_u=zeros(maxIt,1);
order_u=zeros(maxIt,1);
error_v=zeros(maxIt,1);
order_v=zeros(maxIt,1);
error_p=zeros(maxIt,1);
order_p=zeros(maxIt,1);
if ~isfield(pde,'exactphi') || ~isfield(pde,'rhs1')
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

%% Display error and order
% fprintf('\nError_phi\n  dt    &   Error_L2   &  Order \n');
% for k = 1:maxIt
%     fprintf('%.5e  &  %.4e  &  %.2f \n',dt_array(k),error_phi(k),order_phi(k));
% end
% 
% fprintf('\nError_u\n  dt    &   Error_L2   &  Order \n');
% for k = 1:maxIt
%     fprintf('%.5e  &  %.4e  &  %.2f \n',dt_array(k),error_u(k),order_u(k));
% end
% fprintf('\n')
% 
% fprintf('Error_p\n  dt    &   Error_L2   &  Order \n');
% for k = 1:maxIt
%     fprintf('%.5e  &  %.4e  &  %.2f \n',dt_array(k),error_p(k),order_p(k));
% end
% fprintf('\n')
fprintf('   dt \t   & phi_Error_L2 &phi_Order &  u_Error_L2  &  u_Order &  p_Error_L2  &  p_Order \n');
for k = 1:maxIt
    fprintf('%.5e  &  %.4e  &  %.4f  &  %.4e  &  %.4f  &  %.4e  &  %.4f \n',dt_array(k),error_phi(k),order_phi(k),error_u(k),order_u(k),error_p(k),order_p(k));
end
X = ['function:',PDE,'method:',scheme];
disp(X)
fprintf('\n')
%% Plot
figure(2)
hh=loglog(dt_array,error_phi,'*-','LineWidth',3,'MarkerSize',15);
hold on;
hh=loglog(dt_array,error_u, 'p-', 'LineWidth',3,'MarkerSize',20);
hh=loglog(dt_array,error_p, 'b-o','LineWidth',3,'MarkerSize',10);
hh=loglog(dt_array,dt_array,'k-.','LineWidth',3,'MarkerSize',15);

tf = strcmp(scheme2,'_1st');
if tf == 1
    hh=loglog(dt_array,dt_array,'k-.','LineWidth',3,'MarkerSize',10);
    hold on
    h=legend('$L^2$-error: $\phi$','$L^2$-error: $\mathbf{u}$','$L^2$-error: $p$','$\mathcal{O}(\delta t)$','Location','SouthEast');
%     h = legend('Linear scheme: $\phi$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
else 
    hh=loglog(dt_array,dt_array.^2,'k-.','LineWidth',3,'MarkerSize',10);
    hold on
    h = legend('nonLinear scheme: $\phi$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
end
set(h,'interpreter','latex','FontSize',15);
xlabel('Time step $\delta t$','Interpreter','latex');
%  ylabel('Error\_max')
ylabel('$L^2$ error','Interpreter','latex');

ylim([1e-12 5e-1])

grid on;
hold on;

%% Save error and order
name=[scheme,'_phi_e',num2str(para.epsilon),'gamma=',num2str(para.gamma),'S1=',num2str(para.S1),...
      'Nx=',num2str(N),'Ny=',num2str(N)];
% fileID = fopen([name,'.txt'],'w');
% % fprintf(fileID,'%6s\n','%% Results');
% % fprintf(fileID,'%6s\n','% dt	   &   Error_L2	   &  Order');
% % A = [dt_array error];
% % fprintf(fileID,'%.12f   %.4e   \n',A');
% fprintf(fileID,'%.12f     %.4e      %.2f \n',[dt_array,error,order]');
% fclose(fileID);
fileID = fopen([name,'.txt'],'w');
fprintf(fileID,'%.12f   %.4e    %.2f  %.4e   %.2f  %.4e   %.2f\n',...
     [dt_array,error_phi,order_phi,error_u,order_u,error_p,order_p]');
fclose(fileID);

%% results:
