close all;
clear; clc;

% add path
addpath('../','-begin');

%% Space: Domain and N
% domain.left   = -pi;
% domain.right  =  pi;
% domain.bottom = -pi;
% domain.top    =  pi;

domain.xa   = 0;
domain.xb   = 2*pi;
domain.ya   = 0;
domain.yb   = 2*pi;
domain.za   = 0;
domain.zb   = 2*pi;

Lx = domain.xb - domain.xa;
Ly = domain.yb - domain.ya;
Lz = domain.zb - domain.za;

N  = 64;
Nx = N;
Ny = N;
Nz = N;

for tk = 1:2
    if tk == 1
        scheme1 = 'linear';
    else
        scheme1 = 'nonlinear';
    end

    for kt = 1:2
        if kt == 1
            scheme2 = '_1st';   % First-order scheme
%             scheme2 = '_2cn';   % Second-order CN scheme
        else
            scheme2 = '_bdf2';  % Second-order BDF scheme
        end
% scheme1 = 'MSAV';
% scheme2 = '_bdf2';

scheme = [scheme1 , scheme2];


% Parameters
para.epsilon = 6*pi/128;
% para.epsilon = 0.08;
para.gamma = 1;
para.Re = 1;
% para.lambda = 0.001;
para.lambda = 1;
para.S1 = 4;
para.S2 = 4;
para.S3 = 1;

PDE = 'data1';

if 1 == strcmp(PDE,'data1')
    T = 0.02;
    dt_array = 0.001./2.^(1:6)';
    dt_ref = 1e-6;
    para.C0 = 10000; % SAV
    pde = ex03_1_Vesicles_data(para);
end

if 1 == strcmp(PDE,'data2')
    T = 0.002;
    dt_array = 0.001./2.^(6:11)';
    dt_ref = 1e-7;
    para.C0 = 10000;
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
    solver_fun = @CAC_Vesicle_with_NS_3D_linear_newLM_p_SAV_1st;
elseif 1 == strcmp(scheme,'linear_2cn')
    solver_fun = @CAC_Vesicle_3D_linear_newLM_p_SAV_2cn;
elseif 1 == strcmp(scheme,'linear_bdf2')
    solver_fun = @CAC_Vesicle_with_NS_3D_linear_newLM_p_SAV_bdf2;
elseif 1 == strcmp(scheme,'nonlinear_1st')
    solver_fun = @CAC_Vesicle_with_NS_3D_newLM_p_SAV_1st;
elseif 1 == strcmp(scheme,'nonlinear_2cn')
    solver_fun = @CAC_Vesicle_3D_newLM_p_SAV_2cn;
elseif 1 == strcmp(scheme,'nonlinear_bdf2')
    solver_fun = @CAC_Vesicle_with_NS_3D_newLM_p_SAV_bdf2;
elseif 1 == strcmp(scheme,'MSAV_1st')
    solver_fun = @CAC_Vesicle_3D_MSAV_1st;
elseif 1 == strcmp(scheme,'MSAV_bdf2')
    solver_fun = @CAC_Vesicle_3D_MSAV_bdf2;
end


%% Run:
% delete *.mat
if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
    time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
    solver_fun(pde,domain,Nx,Ny,Nz,time,option);
end
for k = 1:maxIt
    dt = dt_array(k);
    time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
    v2 = solver_fun(pde,domain,Nx,Ny,Nz,time,option);
end


%% Compute order of convergence
% error=zeros(maxIt,1);
% order=zeros(maxIt,1);
% if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
%     name=[scheme,'_phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
%           'Nx=',num2str(Nx),'Ny=',num2str(Ny),'Nz=',num2str(Nz),'dt=',num2str(dt_ref)];
%     filename=[name '.mat'];
%     load(filename,'phi');
%     phi_exact = phi;
%     clear phi;
% else
%     Lx = domain.xb - domain.xa;
%     Ly = domain.yb - domain.ya;
%     Lz = domain.zb - domain.za;
%     hx = Lx/Nx;
%     hy = Ly/Ny;
%     hz = Lz/Nz;
%     x  = domain.xa + hx*(0:Nx-1);
%     y  = domain.ya + hy*(0:Ny-1);
%     z  = domain.za + hz*(0:Nz-1);
%     [xx,yy,zz] = ndgrid(x,y,z);
%     phi_exact = pde.exact(xx,yy,zz,T);
% end
% for k = 1:maxIt
%     dt = dt_array(k);
%     name=[scheme,'_phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
%         'Nx=',num2str(Nx),'Ny=',num2str(Ny),'Nz=',num2str(Nz),'dt=',num2str(dt)];
%     filenamek=[name '.mat'];
%     load(filenamek,'phi','hx','hy','hz');
%     err    = fft2((phi_exact - phi).^2);
%     error(k,1) = sqrt(err(1,1,1)*hx*hy);   % L2
%     clear phi;
% end
% order(2:maxIt) = log(error(1:maxIt-1)./error(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));
% 
% %% Display error and order
% fprintf('dt     &   Error_L2\t   &  Order \n');
% for k = 1:maxIt
%     fprintf('%.4e  %.6e  %.2f \n',dt_array(k),error(k),order(k));
%     %         fprintf('1/%d\t& %.6e\t& %.4f %s \n',1./dt_array(k),error(k),order(k),'\\');
% end
% fprintf('\n')
% 
% %% Plot
% figure(2)
% hh=loglog(dt_array,error,'*-','LineWidth',3,'MarkerSize',10);
% hold on 
% tf = strcmp(scheme2,'_1st');
% if tf == 1
%     hh=loglog(dt_array,dt_array,'k-.','LineWidth',3,'MarkerSize',10);
%     hold on
%     h = legend('Linear scheme: $\phi$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
% else 
%     hh=loglog(dt_array,dt_array.^2,'k-.','LineWidth',3,'MarkerSize',10);
%     hold on
%     h = legend('nonLinear scheme: $\phi$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
% end
% set(h,'interpreter','latex','FontSize',15);
% xlabel('Time step $\delta t$','Interpreter','latex');
% %  ylabel('Error\_max')
% ylabel('$L^2$ error','Interpreter','latex');
% grid on;
% hold on;
% 
% %% Save error and order
% name=[scheme,'_phi_e',num2str(para.epsilon),'M',num2str(para.M),...
%       'Nx=',num2str(N),'Ny=',num2str(N),'Nz=',num2str(N)];
% fileID = fopen([name,'.txt'],'w');
% % fprintf(fileID,'%6s\n','%% Results');
% % fprintf(fileID,'%6s\n','% dt	   &   Error_L2	   &  Order');
% % A = [dt_array error];
% % fprintf(fileID,'%.12f   %.4e   \n',A');
% fprintf(fileID,'%.12f     %.4e      %.2f \n',[dt_array,error,order]');
% fclose(fileID);

    end
end

%% results:

%% MSAV
% epsilon=0.147,t=0.0020/0.002, dt=1.56e-06, Nx=64, Ny=64, Nz=64, timeElapsed=151.186315
% dt     &   Error_L2	   &  Order 
% 5.0000e-05  5.449284e-05  0.00 
% 2.5000e-05  1.377393e-05  1.98 
% 1.2500e-05  3.468183e-06  1.99 
% 6.2500e-06  8.703747e-07  1.99 
% 3.1250e-06  2.178878e-07  2.00 
% 1.5625e-06  5.436351e-08  2.00 

%% linear_1st
% ex03_1_Vesicles_data
% lambda=-3.3884e+01,epsilon=0.147,t=0.00200/0.0020, dt=1.56e-06, Nx=64, Ny=64, Nz=64, timeElapsed=442.008188
% dt     &   Error_L2	   &  Order 
% 5.0000e-05  2.646042e-05  0.00 
% 2.5000e-05  1.415568e-05  0.90 
% 1.2500e-05  7.312962e-06  0.95 
% 6.2500e-06  3.696069e-06  0.98 
% 3.1250e-06  1.835461e-06  1.01 
% 1.5625e-06  8.916778e-07  1.04 

%ex03_2_vesicles_data
% dt     &   Error_L2	   &  Order 
% 1.5625e-05  5.390132e-05  0.00 
% 7.8125e-06  2.718519e-05  0.99 
% 3.9063e-06  1.350937e-05  1.01 
% 1.9531e-06  6.597901e-06  1.03 
% 9.7656e-07  3.125411e-06  1.08 
% 4.8828e-07  1.385353e-06  1.17 


%% linear_bdf2
% ex03_1_Vesicles_data
% lambda=-3.3884e+01,epsilon=0.147,t=0.00200/0.0020, dt=1.56e-06, Nx=64, Ny=64, Nz=64,timeElapsed=437.216347
% dt     &   Error_L2	   &  Order 
% 5.0000e-05  1.067610e-06  0.00 
% 2.5000e-05  2.728146e-07  1.97 
% 1.2500e-05  6.903499e-08  1.98 
% 6.2500e-06  1.736633e-08  1.99 
% 3.1250e-06  4.352376e-09  2.00 
% 1.5625e-06  1.086370e-09  2.00 

%% nonlinear_1st
% lambda=-3.3884e+01,epsilon=0.147,t=0.00200/0.0020, dt=1.56e-06, Nx=64, Ny=64, Nz=64, timeElapsed=969.662735
% dt     &   Error_L2	   &  Order 
% 5.0000e-05  2.646042e-05  0.00 
% 2.5000e-05  1.415568e-05  0.90 
% 1.2500e-05  7.312961e-06  0.95 
% 6.2500e-06  3.696069e-06  0.98 
% 3.1250e-06  1.835461e-06  1.01 
% 1.5625e-06  8.916777e-07  1.04 

%% nonlinear_bdf2
% ex03_1_Vesicles_data
% eta=1.0000e+00,lambda=-3.3884e+01,epsilon=0.147,t=0.00200/0.0020, dt=1.56e-06, Nx=64, Ny=64, Nz=64,timeElapsed=737.862854
% dt     &   Error_L2	   &  Order 
% 5.0000e-05  1.067610e-06  0.00 
% 2.5000e-05  2.728146e-07  1.97 
% 1.2500e-05  6.903499e-08  1.98 
% 6.2500e-06  1.736634e-08  1.99 
% 3.1250e-06  4.352379e-09  2.00 
% 1.5625e-06  1.086373e-09  2.00 

