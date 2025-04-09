close all;
clear; clc;

% add path
addpath('../','-begin');

% Space: Domain and N
domain.xa = 0;
domain.xb = 2*pi;
domain.ya = 0;
domain.yb = 2*pi;
domain.za = 0;
domain.zb = 2*pi;

Lx = domain.xb - domain.xa;
Ly = domain.yb - domain.ya;
Lz = domain.zb - domain.za;

% Parameters
para.C0 = 100; % SAV
para.epsilon = 6*pi/128;
para.M = 0.03;

para.name = 'ex19_Vesicles_data';

PDE = 'data2';

if 1 == strcmp(PDE,'data2')
    N  = 64;
    T = 2;
    array = 0:4;
    dt_array = 0.01./2.^array';
    pde = ex02_Vesicles_data(para);
end

Nx = N;
Ny = N;
Nz = N;

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

%% Run:
% delete *.mat
% if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
%     for k = 1:maxIt
%         dt = dt_array(k);
%         time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
%         CAC_Vesicle_3D_BDF2_newLM_p_SAV(pde,domain,Nx,Ny,Nz,time,option);
%     end
% end

%% energy 
figure;
hold on;
lineType = {'-','-','-','-' ,'-','b-','-','-','-','-k','-b'};
for k = 1:maxIt
    energy=load([pde.name,num2str(dt_array(k)),'_energy.txt']);
    tmp = 1;
    plot(energy(tmp:1:end,1),energy(tmp:1:end,2),char(lineType(k)),'LineWidth',2.5);
    legend_str{k} = ['$',num2str(k),':\delta t = 0.01/2^',num2str(array(k)),'$'];
end
h = legend(legend_str,'box','off','interpreter','latex');
% set(h,'interpreter','latex');
xlabel('Time','Fontsize',24);ylabel('Original Energy','Fontsize',24,'interpreter','latex');
% xlim([0,5])
% % ylim([6490,6650])
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

figname1 = ['C:\Users\heyan\Desktop\sjj\songfig\Original_energy_3dvesicle','.png'];
print(figname1,'-dpng', '-r300')



%% mass and surface
figure;
hold on;
% lineType = {'.-', 's-','*-','o-','+-' ,'.-' ,'--k',':','-b'};
% lineType = {'>-', 's-','*-','o-','+-' ,'.-' ,'--k','-b'};
% lineType = {'>-', 's-','*-','o-','+-' ,'.-' ,'--k','-b'};

lineType = {'-','-','-','-' ,'-','b-','-','-','-','-k','-b'};
for k = 1:maxIt
    mass=load([pde.name,num2str(dt_array(k)),'_mass.txt']);
%     energy=load([pde.name,'_dt_',num2str(dt_array(k)),'_mass.txt']);
    tmp = 1;
    plot(mass(tmp:1:end,1),abs((mass(tmp:1:end,3)-mass(1,3))./mass(1,3)),char(lineType(k)),'LineWidth',2.5);
    legend_str{k} = ['$',num2str(k),':\delta t = 0.01/2^',num2str(array(k)),'$'];
end
h = legend(legend_str,'box','off','interpreter','latex');
% set(h,'interpreter','latex');
xlabel('Time','Fontsize',24);ylabel('Ratio of Surface Area Difference','Fontsize',0.5,'interpreter','latex');
% xlim([0,5])
% % ylim([6490,6650]
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

figname2 = ['C:\Users\heyan\Desktop\sjj\songfig\surface_ratio_3dvesicle','.png'];
print(figname2,'-dpng', '-r300')

