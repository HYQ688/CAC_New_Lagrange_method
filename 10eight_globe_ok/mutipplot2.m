close all;
clear; clc;

% add path
addpath('../','-begin');

% Space: Domain and N
domain.xa =   0 ;
domain.xb =2*pi ;
domain.ya =   0 ;
domain.yb =2*pi ;
domain.za =   0 ;
domain.zb =2*pi ;

Lx = domain.xb - domain.xa;
Ly = domain.yb - domain.ya;
Lz = domain.zb - domain.za;

N = 64;
Nx = N;
Ny = N;
Nz = N;

hx = Lx/Nx;
hy = Ly/Ny;
hz = Lz/Nz;
x  = domain.xa + hx*(0:Nx-1);
y  = domain.ya + hy*(0:Ny-1);
z  = domain.za + hz*(0:Nz-1);

[xx,yy,zz] = ndgrid(x,y,z);

% Parameters
para.epsilon = 6*pi/128;
para.gamma = 1;
para.Re = 1;
para.lambda = 0.1;
para.C0 = 10000;
para.S1 = 4;
para.S2 = 4;
para.S3 = 1;

%     sch1 = 'linear';
%     sch1 = 'nonlinear';
    sch1 = 'MSAV';

%     sch2 = '_1st'; 
    sch2 = '_bdf2';

    sch3 = 'linear';
%     sch3 = 'nonlinear';
%     sch3 = 'MSAV';

%     sch4 = '_1st'; 
    sch4 = '_bdf2';

    sch5 = 'nonlinear';
    sch6 = '_bdf2';
    scheme1 = [sch1 , sch2];
    scheme2 = [sch3 , sch4];
    scheme3 = [sch5,sch6];

pdename1 = [scheme1,'_ex03_Vesicles_data_dt_';];
pdename2 = [scheme2,'_ex02_Vesicles_data_dt';];
pdename3 = [scheme3,'_ex02_Vesicles_data_dt';];
pde = ex02_Vesicles_data(para);

T =20;
t0 = 0;
tsave = 0.02*T;

dt_array = 0.01./2.^(0:4); 
dt_ref = 0.01/2^3;

maxIt = length(dt_array);

% 能量
lineType = {'b-', 'r--', 'r:'}; 

figure; % 创建新的图形窗口
hold on; % 保持当前绘图，以便添加其他元素
mass1=load([pdename2,num2str(dt_ref),'_mass_data.txt']);
mass2=load([pdename3,num2str(dt_ref),'_mass_data.txt']);
tmp = 1;

fillcolor1=[0.85, 0.33, 0.10];
fillcolor2=[0.93, 0.69, 0.13];
fillcolor3=[0.00, 0.45, 0.74];

    % 绘制第一条线
    h1 = plot(mass1(tmp:1:end, 1), abs((mass1(tmp:1:end, 3)-mass1(1,3))/mass1(1,3)), char(lineType(1)), 'LineWidth', 2.5,'Color',fillcolor3);
    % 绘制第二条线
    h2 = plot(mass2(tmp:1:end, 1), abs((mass2(tmp:1:end, 3)-mass2(1,3))/mass2(1,3)), char(lineType(2)), 'LineWidth', 2.5,'Color',fillcolor1);
%     % 绘制第三条线
%     % h3 = plot(energy(tmp:1:end, 1), energy(tmp:1:end, 4), char(lineType(3)), 'LineWidth', 2.5);
% 
ylim([0,1.5e-4]);
h = legend('LM-SAV','NLM','box','off','interpreter','latex');
set(h,'interpreter','latex','FontSize',15);
xlabel('Time','Fontsize',20);ylabel('Surface','Fontsize',20,'interpreter','latex');

% 绘制能量曲线
 hold on;
box on; % 添加边框
grid on; % 显示网格
set(gca, 'LineWidth', 1.5); % 设置边框线条粗细

% 下面部分创建主图
dirname = '_ex20_BDF2_3D_CAC_Vesicles_data_LM';
dirname = [scheme3,dirname];
datadir = [dirname,'/data'];
X = load([datadir, '/X.txt']);
Y = load([datadir, '/Y.txt']);
Z = load([datadir, '/Z.txt']);

%小图1
t=0;
ssp = [datadir '/phi_t='  num2str(t) '.txt'];
phi = load(ssp);


% 在主图中创建小图
pos1 = [0.01, 0.4, 0.4, 0.4]; % 位置和大小 [x, y, width, height]
ax = axes('Position', pos1); % 创建小图区域
% 使用 axes 创建小图
% 修正：生成三维网格数据

% 确保 phi 是三维数组并且与 x, y, z 的网格尺寸匹配
% 绘制等值面
% p = patch(isosurface(x, y, z, phi, 0.5)); % 绘制等值面
% set(p, 'FaceColor', 'red', 'EdgeColor', 'none'); % 设置颜色
% view(3); % 设置为 3D 视角
% grid on;
% box on;
% title('3D Example', 'FontSize', 12);

    phi = reshape(phi,Nx,Ny,Nz);    
    p0 = patch(isosurface(xx,yy,zz,phi,0));
    set(p0,'FaceColor','red', 'EdgeColor','none');
    % view(3);
    daspect([1 1 1]); 
    camlight;  
    lighting phong;
     % axis square;
    axis off;
    %  xlim([0, 2*pi])
    % ylim([0, 2*pi])
    % zlim([0, 2*pi])
    % box on; 
     view(-161,16);
%     view(115,20)    

    % ax = gca;
    % c = ax.Color;
    % ax.Color = [0.30,0.75,0.90];
     drawnow;

    set(gcf, 'InvertHardCopy', 'off');
    set(0,'defaultfigurecolor','w') 
% axes('Position', pos1);
% s = pcolor(X, Y, phi);
% s.FaceColor = 'interp';
% s.EdgeColor = 'interp';
% colormap jet;
% axis square;
% axis tight;
% axis off;
title(['T=',num2str(t)]);
hold on;

%小图2
t=2;
ssp = [datadir '/phi_t='  num2str(t) '.txt'];
phi = load(ssp);


% 在主图中创建小图
pos1 = [0.18, 0.45, 0.4, 0.4]; % 位置和大小 [x, y, width, height]
ax = axes('Position', pos1); % 创建小图区域
% 使用 axes 创建小图
% 修正：生成三维网格数据

% 确保 phi 是三维数组并且与 x, y, z 的网格尺寸匹配
% 绘制等值面
% p = patch(isosurface(x, y, z, phi, 0.5)); % 绘制等值面
% set(p, 'FaceColor', 'red', 'EdgeColor', 'none'); % 设置颜色
% view(3); % 设置为 3D 视角
% grid on;
% box on;
% title('3D Example', 'FontSize', 12);

    phi = reshape(phi,Nx,Ny,Nz);    
    p0 = patch(isosurface(xx,yy,zz,phi,0));
    set(p0,'FaceColor','red', 'EdgeColor','none');
    % view(3);
    daspect([1 1 1]); 
    camlight;  
    lighting phong;
     % axis square;
    axis off;
    %  xlim([0, 2*pi])
    % ylim([0, 2*pi])
    % zlim([0, 2*pi])
    % box on; 
     view(-161,16);
%     view(115,20)    

    % ax = gca;
    % c = ax.Color;
    % ax.Color = [0.30,0.75,0.90];
     drawnow;

    set(gcf, 'InvertHardCopy', 'off');
    set(0,'defaultfigurecolor','w') 
% axes('Position', pos1);
% s = pcolor(X, Y, phi);
% s.FaceColor = 'interp';
% s.EdgeColor = 'interp';
% colormap jet;
% axis square;
% axis tight;
% axis off;
title(['T=',num2str(t)]);
hold on;

%小图3
t=10;
ssp = [datadir '/phi_t='  num2str(t) '.txt'];
phi = load(ssp);


% 在主图中创建小图
pos1 = [0.38, 0.42, 0.4, 0.4]; % 位置和大小 [x, y, width, height]
ax = axes('Position', pos1); % 创建小图区域
% 使用 axes 创建小图
% 修正：生成三维网格数据

% 确保 phi 是三维数组并且与 x, y, z 的网格尺寸匹配
% 绘制等值面
% p = patch(isosurface(x, y, z, phi, 0.5)); % 绘制等值面
% set(p, 'FaceColor', 'red', 'EdgeColor', 'none'); % 设置颜色
% view(3); % 设置为 3D 视角
% grid on;
% box on;
% title('3D Example', 'FontSize', 12);

    phi = reshape(phi,Nx,Ny,Nz);    
    p0 = patch(isosurface(xx,yy,zz,phi,0));
    set(p0,'FaceColor','red', 'EdgeColor','none');
    % view(3);
    daspect([1 1 1]); 
    camlight;  
    lighting phong;
     % axis square;
    axis off;
    % axis equal;
    %  xlim([0, 2*pi])
    % ylim([0, 2*pi])
    % zlim([0, 2*pi])
    % box on; 
     view(-161,16);
%     view(115,20)    

    % ax = gca;
    % c = ax.Color;
    % ax.Color = [0.30,0.75,0.90];
     drawnow;

    set(gcf, 'InvertHardCopy', 'off');
    set(0,'defaultfigurecolor','w') 
% axes('Position', pos1);
% s = pcolor(X, Y, phi);
% s.FaceColor = 'interp';
% s.EdgeColor = 'interp';
% colormap jet;
% axis square;
% axis tight;
% axis off;
title(['T=',num2str(t)]);
hold on;


%小图4
t=20;
ssp = [datadir '/phi_t='  num2str(t) '.txt'];
phi = load(ssp);


% 在主图中创建小图
pos1 = [0.6, 0.35, 0.4, 0.4]; % 位置和大小 [x, y, width, height]
ax = axes('Position', pos1); % 创建小图区域
% 使用 axes 创建小图
% 修正：生成三维网格数据

% 确保 phi 是三维数组并且与 x, y, z 的网格尺寸匹配
% 绘制等值面
% p = patch(isosurface(x, y, z, phi, 0.5)); % 绘制等值面
% set(p, 'FaceColor', 'red', 'EdgeColor', 'none'); % 设置颜色
% view(3); % 设置为 3D 视角
% grid on;
% box on;
% title('3D Example', 'FontSize', 12);

    phi = reshape(phi,Nx,Ny,Nz);    
    p0 = patch(isosurface(xx,yy,zz,phi,0));
    set(p0,'FaceColor','red', 'EdgeColor','none');
    % view(3);
    daspect([1 1 1]); 
    camlight;  
    lighting phong;
     % axis square;
    axis off;
    % axis equal;
    %  xlim([0, 2*pi])
    % ylim([0, 2*pi])
    % zlim([0, 2*pi])
    % box on; 
     view(-161,16);
%     view(115,20)    

    % ax = gca;
    % c = ax.Color;
    % ax.Color = [0.30,0.75,0.90];
     drawnow;

    set(gcf, 'InvertHardCopy', 'off');
    set(0,'defaultfigurecolor','w') 
% axes('Position', pos1);
% s = pcolor(X, Y, phi);
% s.FaceColor = 'interp';
% s.EdgeColor = 'interp';
% colormap jet;
% axis square;
% axis tight;
% axis off;
title(['T=',num2str(t)]);
hold on;






% 
figname = ['D:\paper\phase-field\CAC_Vesicles_surface_LagrangeMultiplier_SAV\fig\',scheme2,'N4mass.png'];
print(figname,'-dpng', '-r300')








