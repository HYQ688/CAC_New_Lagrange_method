close all;
clear; clc;
tic
% add path
addpath('../','-begin');

% Space: Domain and N
domain.left   = 0;
domain.right  = 2*pi;
domain.bottom = 0;
domain.top    = 2*pi;

Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;

% scheme1 = 'linear';
scheme1 = 'nonlinear';
% scheme1 = 'MSAV';

% scheme2 = '_1st';   % First-order scheme
% scheme2 = '_2cn';   % Second-order CN scheme
scheme2 = '_bdf2';  % Second-order BDF scheme

scheme = [scheme1 , scheme2];


if  1 == strcmp(scheme1,'MSAV')
% MSAV Parameters

para.C0 = 100; % SAV
para.Lx = Lx;
para.Ly = Ly;
para.epsilon = 6*pi/128;
para.M  = 1;
para.M1 = 50;
para.M2 = 100;
para.S1 = 4;
para.S2 = 4;
para.S3 = 1;
para.beta_m = 1;

para.name = [scheme,'_ex03_Vesicles_data'];
pde = ex03_Vesicles_data(para);

else
% Parameters
%linear
% para.epsilon = 6*pi/128;
% para.gamma = 1;
% para.Re = 1;
% para.lambda = 1;
% para.C0 = 100;
% para.S1 = 4;
% para.S2 = 4;
% para.S3 = 1;

%nonlinear
para.epsilon = 6*pi/128;
para.gamma = 1;
para.Re = 1;
para.lambda = 1;
para.C0 = 100;
para.S1 = 4;
para.S2 = 4;
para.S3 = 1;

para.name = [scheme,'_ex02_Vesicles_data'];
pde = ex02_Vesicles_data(para);

end

N  = 64;
% T = 2; %linear_1st  linear_bdf2
T =4; % nonlinear_1st nonlinear_bdf2
% T = 2;
% array = 2:6;
array = 3;
dt_array = 0.01./2.^array';
dt_ref = 1e-5;

Nx = N;
Ny = N;

% Time: dt T
t0 = 0;
tsave = 0.5*T;

maxIt = length(dt_array);

%% option
option.plotflag   = 0;
option.printflag  = 1;
option.vtkflag    = 0;
option.saveflag   = 0;
option.savefinal  = 0;
option.energyflag = 1;
option.tol = 1e-14;
option.tolit = 1e-11;
option.maxit = 2000;

if 1 == strcmp(scheme,'linear_1st')
    solver_fun = @CAC_Vesicle_with_NS_2D_linear_newLM_p_SAV_1st;
elseif 1 == strcmp(scheme,'linear_2cn')
    solver_fun = @CAC_Vesicle_with_NS_2D_linear_newLM_p_SAV_bdf2;
elseif 1 == strcmp(scheme,'linear_bdf2')
    solver_fun = @CAC_Vesicle_with_NS_2D_linear_newLM_p_SAV_bdf2;
elseif 1 == strcmp(scheme,'nonlinear_1st')
    solver_fun = @CAC_Vesicle_with_NS_2D_newLM_p_SAV_1st;
elseif 1 == strcmp(scheme,'nonlinear_2cn')
    solver_fun = @CAC_Vesicle_2D_newLM_p_SAV_2cn;
elseif 1 == strcmp(scheme,'nonlinear_bdf2')
    solver_fun = @CAC_Vesicle_with_NS_2D_newLM_p_SAV_bdf2;
elseif 1 == strcmp(scheme,'MSAV_1st')
    solver_fun = @CAC_Vesicle_2D_MSAV_1st;
elseif 1 == strcmp(scheme,'MSAV_bdf2')
    solver_fun = @CAC_Vesicle_2D_MSAV_bdf2;
end

%% Run:
% delete *.mat
% if ~isfield(pde,'exactu') || ~isfield(pde,'rhs2')
%     for k = 1:maxIt
%         dt = dt_array(k);
%         time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
%         solver_fun(pde,domain,Nx,Ny,time,option);
%     end
% end
toc
%% Original energy
figure(3);
hold on;
lineType = {'-','-','-','-' ,'-','b-','-','-','-','-k','-b'};
for k = 1:maxIt
    energy=load([pde.name, '_dt',num2str(dt_array(k)),'_S1',num2str(pde.S1), '_energy.txt']);
    tmp = 1;
    plot(energy(tmp:1:end,1),energy(tmp:1:end,2),char(lineType(k)),'LineWidth',3.5);
    legend_str{k} = ['$',num2str(k),':\delta t = 0.01/2^',num2str(array(k)),'$'];
end
h = legend(legend_str,'box','off','interpreter','latex');
set(h,'interpreter','latex');
xlabel('Time','Fontsize',24);ylabel('Original Energy','Fontsize',24,'interpreter','latex');
% xlim([0,5])
ylim([3,16])
% % yticks([0:2:20])
set(gca,'FontSize',22);
set(gca,'linewidth',1.8)
% set(gca,'xtick',0:0.2:T);
% set(gca,'ytick',3:0.1:4);
grid on;
box on;
figure_FontSize=24;
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',figure_FontSize);
set(findobj('FontSize',10),'FontSize',figure_FontSize); %这4句是将字体大小改为8号字，在小图里很清晰

figname1 = ['D:\paper\phase-field\CAC_Vesicles_surface_LagrangeMultiplier_SAV\fig\',scheme,'_2D_Original energy','.png'];
% print(figname1,'-dpng', '-r300')

%% Modified energy
% figure(3);
% hold on;
% lineType = {'-','-','-','-' ,'-','b-','-','-','-','-k','-b'};
% for k = 1:maxIt
%     energy=load([pde.name, '_dt',num2str(dt_array(k)),'_S1',num2str(pde.S1), '_energy.txt']);
%     tmp = 1;
%     plot(energy(tmp:1:end,1),energy(tmp:1:end,3),char(lineType(k)),'LineWidth',2.5);
%     legend_str{k} = ['$',num2str(k),':\delta t = 0.01/2^',num2str(array(k)),'$'];
% end
% h = legend(legend_str,'box','off','interpreter','latex');
% set(h,'interpreter','latex');
% % legend_str = ["Original energy", "Discrete modified energy"];
% % h = legend(legend_str,'box','off','interpreter','latex');
% xlabel('Time','Fontsize',24);ylabel('Discrete modified energy','Fontsize',12,'interpreter','latex');
% % xlim([0,5])
% ylim([3,16])
% % % yticks([0:2:20])
% set(gca,'FontSize',22);
% set(gca,'linewidth',1.8)
% % set(gca,'xtick',0:0.2:T);
% % set(gca,'ytick',3:0.1:4);
% grid on;
% box on;
% figure_FontSize=24;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize);
% set(get(gca,'YLabel'),'FontSize',figure_FontSize);
% set(findobj('FontSize',10),'FontSize',figure_FontSize); %这4句是将字体大小改为8号字，在小图里很清晰
% 
% figname1 = ['D:\paper\phase-field\CAC_Vesicles_surface_LagrangeMultiplier_SAV\fig\',scheme,'_2D_Menergy','.png'];
% % print(figname1,'-dpng', '-r300')

%% discrete modified energy
figure(3);
hold on;
lineType = {'--','-','-','-' ,'-','b-','-','-','-','-k','-b'};
for k = 1:maxIt
    energy=load([pde.name, '_dt',num2str(dt_array(k)),'_S1',num2str(pde.S1), '_energy.txt']);
    tmp = 1;
    plot(energy(tmp:1:end,1),energy(tmp:1:end,3),char(lineType(k)),'LineWidth',2.5);
end
legend_str = ["Original energy",  "Discrete modified energy"];
h = legend(legend_str,'box','off','interpreter','latex');
% set(h,'interpreter','latex');

xlabel('Time','Fontsize',24);ylabel('Energy','Fontsize',24,'interpreter','latex');
% xlim([0,5])
ylim([3,16])
% % yticks([0:2:20])
set(gca,'FontSize',22);
set(gca,'linewidth',1.8)
% set(gca,'xtick',0:0.2:T);
% set(gca,'ytick',3:0.1:4);
grid on;
box on;
figure_FontSize=24;
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',figure_FontSize);
set(findobj('FontSize',10),'FontSize',figure_FontSize); %这4句是将字体大小改为8号字，在小图里很清晰

figname1 = ['D:\paper\phase-field\CAC_Vesicles_surface_LagrangeMultiplier_SAV\fig\',scheme,'_2D_Cenergy','.png'];
% print(figname1,'-dpng', '-r300')

%% mass and surface
figure;
hold on;
% lineType = {'.-', 's-','*-','o-','+-' ,'.-' ,'--k',':','-b'};
% lineType = {'>-', 's-','*-','o-','+-' ,'.-' ,'--k','-b'};
% lineType = {'>-', 's-','*-','o-','+-' ,'.-' ,'--k','-b'};

lineType = {'-','-','-','-' ,'-','b-','-','-','-','-k','-b'};
for k = 1:maxIt
    mass=load([pde.name,'_dt',num2str(dt_array(k)),'_S1',num2str(pde.S1),'_mass_data.txt']);
%     energy=load([pde.name,'_dt_',num2str(dt_array(k)),'_mass.txt']);
    tmp = 1;
    plot(mass(tmp:1:end,1),abs((mass(tmp:1:end,3)-mass(1,3))./mass(1,3)),char(lineType(k)),'LineWidth',2.5);
    legend_str{k} = ['$',num2str(k),':\delta t = 0.01/2^',num2str(array(k)),'$'];
end
h = legend(legend_str,'box','off','interpreter','latex');
% set(h,'interpreter','latex');
xlabel('Time','Fontsize',24);ylabel('Ratio of Surface Area Difference','Fontsize',0.5,'interpreter','latex');
% xlim([0,5])
ylim([0,4e-4])
% % yticks([0:2:20])
set(gca,'FontSize',22);
set(gca,'linewidth',1.8)
% set(gca,'xtick',0:0.2:T);
% set(gca,'ytick',3:0.1:4);
grid on;
box on;

figure_FontSize=24;
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',18);
set(findobj('FontSize',10),'FontSize',figure_FontSize); %这4句是将字体大小改为8号字，在小图里很清晰

figname2 = ['D:\paper\phase-field\CAC_Vesicles_surface_LagrangeMultiplier_SAV\fig\',scheme,'_2D_surface','.png'];
% print(figname2,'-dpng', '-r300')





