close all;
clear; clc;

% add path
addpath('../','-begin');

% Space: Domain and N
domain.left   =   0 ;
domain.right  =2*pi ;
domain.bottom =   0 ;
domain.top    =2*pi ;

Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;

N = 128;
Nx = N;
Ny = N;

% Parameters
para.epsilon = 6*pi/128;
para.gamma = 1;
para.Re = 1;
para.lambda = 0.1;
para.C0 = 10000;
para.S1 = 4;
para.S2 = 4;
para.S3 = 1;

% scheme1 = 'linear';
scheme1 = 'nonlinear';
scheme2 = '_bdf2';  % Second-order BDF scheme
scheme = [scheme1 , scheme2];

para.name = [scheme,'_ex02_Vesicles_data'];
pde = ex02_Vesicles_data(para);

T = 4;
t0 = 0;
tsave = 0.02*T;

dt_array = 0.01./2.^(0:4); 
dt_ref = 0.01/2^3;

maxIt = length(dt_array);

% 能量
lineType = {'b-', 'g--', 'r:'}; 

figure; % 创建新的图形窗口
hold on; % 保持当前绘图，以便添加其他元素

for k = 1:maxIt
    energy = load([pde.name, '_dt0.00125_S14_energy.txt']);
    tmp = 1;

    % 绘制第一条线
    h1 = plot(energy(tmp:1:end, 1), energy(tmp:1:end, 2), char(lineType(1)), 'LineWidth', 2.5);
    % 绘制第二条线
    h2 = plot(energy(tmp:1:end, 1), energy(tmp:1:end, 3), char(lineType(2)), 'LineWidth', 2.5);
    % 绘制第三条线
    % h3 = plot(energy(tmp:1:end, 1), energy(tmp:1:end, 4), char(lineType(3)), 'LineWidth', 2.5);
end

% 定义图例字符串
% legend_str = {'Original Energy', 'Modified Energy', 'Discrete modified Energy'};
% h = legend([h1, h2, h3], legend_str, 'box', 'off', 'Interpreter', 'latex', 'Location', 'north east');

legend_str = {'Original Energy',  'Discrete Energy'};
h = legend([h1, h2], legend_str, 'box', 'off', 'Interpreter', 'latex', 'Location', 'north east');

% 设置标签字体和大小

xlabel('Time','Fontsize',24);
ylabel('Energy','Fontsize',24,'interpreter','latex');
ylim([0.2,2.1]);
set(gca,'FontSize',12);
grid on;
box on;
set(gca,'linewidth',2.5);

% 下面部分创建主图
dirname = '_CAC_Vesicles_data_LM_six_circles';
dirname = [scheme,dirname];
datadir = [dirname,'/data'];
X = load([datadir, '/X.txt']);
Y = load([datadir, '/Y.txt']);

%小图1
t=0;
ssp = [datadir '/phi_t='  num2str(t) '.txt'];
phi = load(ssp);
    
% 在主图中创建小图
pos1 = [0.2, 0.62, 0.24, 0.24]; % 位置和大小 [x, y, width, height]

% 使用 axes 创建小图
axes('Position', pos1);
s = pcolor(X, Y, phi);
s.FaceColor = 'interp';
s.EdgeColor = 'interp';
colormap jet;
axis square;
axis tight;
axis off;
title(['T=',num2str(t)]);
hold on;

% 设置其他图形参数
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1.8);

annotation('arrow', [0.24, 0.13], [0.73, 0.9]);

%小图2
t=2;
ssp = [datadir '/phi_t='  num2str(t) '.txt'];
phi = load(ssp);
    
% 在主图中创建小图
pos1 = [0.25, 0.31, 0.24, 0.24]; % 位置和大小 [x, y, width, height]

% 使用 axes 创建小图
axes('Position', pos1);
s = pcolor(X, Y, phi);
s.FaceColor = 'interp';
s.EdgeColor = 'interp';
colormap jet;
axis square;
axis tight;
axis off;
title(['T=',num2str(t)]);
hold on;

% 设置其他图形参数
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1.8);

annotation('arrow', [0.28, 0.21], [0.42, 0.365]);

%小图3
t=10;
ssp = [datadir '/phi_t='  num2str(t) '.txt'];
phi = load(ssp);
    
% 在主图中创建小图
pos1 = [0.47, 0.28, 0.25, 0.25]; % 位置和大小 [x, y, width, height]

% 使用 axes 创建小图
axes('Position', pos1);
s = pcolor(X, Y, phi);
s.FaceColor = 'interp';
s.EdgeColor = 'interp';
colormap jet;
axis square;
axis tight;
axis off;
title(['T=',num2str(t)]);
hold on;

% 设置其他图形参数
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1.8);

annotation('arrow', [0.58, 0.52], [0.28, 0.24]);

%小图4
t=20;
ssp = [datadir '/phi_t='  num2str(t) '.txt'];
phi = load(ssp);
    
% 在主图中创建小图
pos1 = [0.68, 0.28, 0.25, 0.25]; % 位置和大小 [x, y, width, height]

% 使用 axes 创建小图
axes('Position', pos1);
s = pcolor(X, Y, phi);
s.FaceColor = 'interp';
s.EdgeColor = 'interp';
colormap jet;
axis square;
axis tight;
axis off;
title(['T=',num2str(t)]);
hold on;

% 设置其他图形参数
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1.8);

annotation('arrow', [0.81, 0.9], [0.28, 0.2]);

figname = ['D:\paper\phase-field\CAC_Vesicles_surface_LagrangeMultiplier_SAV\fig\',scheme,'N1energy.png'];
% print(figname,'-dpng', '-r300')









