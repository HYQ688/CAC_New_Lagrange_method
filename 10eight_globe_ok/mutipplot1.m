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
energy1=load([pdename2,num2str(dt_ref),'_energy.txt']);
energy2=load([pdename3,num2str(dt_ref),'_energy.txt']);
tmp = 1;

fillcolor1=[0.85, 0.33, 0.10];
fillcolor2=[0.93, 0.69, 0.13];
fillcolor3=[0.00, 0.45, 0.74];

    % 绘制第一条线
    h1 = plot(energy1(tmp:1:end, 1), energy1(tmp:1:end, 2), char(lineType(1)), 'LineWidth', 2.5,'Color',fillcolor3);
    % 绘制第二条线
    h2 = plot(energy2(tmp:1:end, 1), energy2(tmp:1:end, 2), char(lineType(2)), 'LineWidth', 2.5,'Color',fillcolor1);
%     % 绘制第三条线
%     % h3 = plot(energy(tmp:1:end, 1), energy(tmp:1:end, 4), char(lineType(3)), 'LineWidth', 2.5);
% 

h = legend('LM-SAV','NLM','box','off','interpreter','latex');
set(h,'interpreter','latex','FontSize',15);
xlabel('Time','Fontsize',20);ylabel('Original Energy','Fontsize',20,'interpreter','latex');

% 绘制能量曲线
 hold on;
box on; % 添加边框
grid on; % 显示网格
set(gca, 'LineWidth', 1.5); % 设置边框线条粗细

% 下面部分创建主图
dirname = '_ex20_BDF2_3D_CAC_Vesicles_data_LM';
dirname = [scheme2,dirname];
datadir = [dirname,'/data'];
X = load([datadir, '/X.txt']);
Y = load([datadir, '/Y.txt']);
Z = load([datadir, '/Z.txt']);

%小图1
t=0;
ssp = [datadir '/phi_t='  num2str(t) '.txt'];
phi = load(ssp);


% 在主图中创建小图
pos1 = [0.05, 0.53, 0.4, 0.4]; % 位置和大小 [x, y, width, height]
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
t = title(['T=', num2str(t)]);
        t.Position = [8, -15, 0]; % 调整标题的位置
hold on;

%小图2
t=2;
ssp = [datadir '/phi_t='  num2str(t) '.txt'];
phi = load(ssp);


% 在主图中创建小图
pos1 = [0.14, 0.25, 0.4, 0.4]; % 位置和大小 [x, y, width, height]
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
t = title(['T=', num2str(t)]);
        t.Position = [8, -15, 0]; % 调整标题的位置
hold on;

%小图3
t=10;
ssp = [datadir '/phi_t='  num2str(t) '.txt'];
phi = load(ssp);


% 在主图中创建小图
pos1 = [0.35, 0.15, 0.4, 0.4]; % 位置和大小 [x, y, width, height]
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
t = title(['T=', num2str(t)]);
        t.Position = [8, -15, 0]; % 调整标题的位置
hold on;


%小图4
t=20;
ssp = [datadir '/phi_t='  num2str(t) '.txt'];
phi = load(ssp);


% 在主图中创建小图
pos1 = [0.58, 0.12, 0.4, 0.4]; % 位置和大小 [x, y, width, height]
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
t = title(['T=', num2str(t)]);
        t.Position = [8, -15, 0]; % 调整标题的位置
hold on;

% 
% 
% 
% 
% 
% 
% close all;
% clear; clc;
% 
% % 添加路径
% addpath('../','-begin');
% 
% % 空间域定义
% domain = struct('xa', 0, 'xb', 2*pi, 'ya', 0, 'yb', 2*pi, 'za', 0, 'zb', 2*pi);
% Lx = domain.xb - domain.xa;
% Ly = domain.yb - domain.ya;
% Lz = domain.zb - domain.za;
% 
% % 网格划分
% N = 64;
% hx = Lx/N; hy = Ly/N; hz = Lz/N;
% x = domain.xa + hx*(0:N-1);
% y = domain.ya + hy*(0:N-1);
% z = domain.za + hz*(0:N-1);
% [xx, yy, zz] = ndgrid(x, y, z);
% 
% % 物理参数
% para = struct('epsilon', 6*pi/128, 'gamma', 1, 'Re', 1, ...
%               'lambda', 0.1, 'C0', 10000, 'S1', 4, 'S2', 4, 'S3', 1);
% 
% % 方案定义
% scheme1 = 'MSAV_bdf2';
% scheme2 = 'linear_bdf2';
% scheme3 = 'nonlinear_bdf2';
% 
% % PDE 文件名
% pdename2 = [scheme2, '_ex02_Vesicles_data_dt'];
% pdename3 = [scheme3, '_ex02_Vesicles_data_dt'];
% 
% % 计算时间步长
% dt_ref = 0.01/2^3;
% 
% % 读取能量数据
% energy1 = load([pdename2, num2str(dt_ref), '_S14_energy.txt']);
% energy2 = load([pdename3, num2str(dt_ref), '_S14_energy.txt']);
% 
% % 绘制能量曲线
% figure; hold on;
% box on; % 添加边框
% grid on; % 显示网格
% set(gca, 'LineWidth', 1.5); % 设置边框线条粗细
% 
% fillcolor1 = [0.85, 0.33, 0.10]; % 颜色定义
% fillcolor3 = [0.00, 0.45, 0.74];
% 
% plot(energy1(:,1), energy1(:,2), 'b-', 'LineWidth', 2.5, 'Color', fillcolor3);
% plot(energy2(:,1), energy2(:,2), 'r--', 'LineWidth', 2.5, 'Color', fillcolor1);
% 
% legend({'LM-SAV', 'NLM'}, 'Interpreter', 'latex', 'FontSize', 15);
% xlabel('Time', 'FontSize', 20);
% ylabel('Original Energy', 'FontSize', 20, 'Interpreter', 'latex');
% 
% % 数据目录
% 
% datadir = [scheme2, '_ex20_BDF2_3D_CAC_Vesicles_data_LM', '/data'];
% X = load(fullfile(datadir, 'X.txt'));
% Y = load(fullfile(datadir, 'Y.txt'));
% Z = load(fullfile(datadir, 'Z.txt'));
% 
% % 需要绘制等值面的时间点
% time_list = [0, 2, 10, 20];
% pos_list = [0.06, 0.54; 0.16, 0.21; 0.38, 0.15; 0.6, 0.12]; % 子图位置
% 
% % 统一绘制 3D 等值面
% for i = 1:length(time_list)
%     t = time_list(i);
%     filename = fullfile(datadir, ['phi_t=', num2str(t), '.txt']);
% 
%     if exist(filename, 'file')
%         phi = load(filename);
%         phi = reshape(phi, N, N, N);
% 
%         % 创建小图
%         axes('Position', [pos_list(i,:), 0.4, 0.4]); 
%         p = patch(isosurface(xx, yy, zz, phi, 0));
%         set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
%         daspect([1,1,1]); camlight; lighting phong; axis off;
%         view(-161, 16);
%         t = title(['T=', num2str(t)]);
%         t.Position = [8, -15, 0]; % 调整标题的位置
%     else
%         warning('File not found: %s', filename);
%     end
% end
% 
% % 设置背景颜色
% set(gcf, 'InvertHardCopy', 'off');
% set(0, 'defaultfigurecolor', 'w');


figname = ['D:\paper\phase-field\CAC_Vesicles_surface_LagrangeMultiplier_SAV\fig\',scheme2,'N4energy.png'];
print(figname,'-dpng', '-r300')

% 
% 
% 
% 
% 
% 
% 
