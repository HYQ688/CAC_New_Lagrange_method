% close all;
% clear; clc;
% 
% % add path
% addpath('../','-begin');
% 
% % Space: Domain and N
% domain.left   =   0 ;
% domain.right  =2*pi ;
% domain.bottom =   0 ;
% domain.top    =2*pi ;
% 
% Lx = domain.right - domain.left;
% Ly = domain.top   - domain.bottom;
% 
% N = 128;
% Nx = N;
% Ny = N;
% 
% % Parameters
% para.epsilon = 6*pi/128;
% para.gamma = 0.1;
% para.Re = 1;
% para.lambda = 0.1;
% para.C0 = 10000;
% para.S1 = 4;
% para.S2 = 4;
% para.S3 = 1;
% 
% % scheme1 = 'linear';
% scheme1 = 'nonlinear';
% scheme2 = '_bdf2';  % Second-order BDF scheme
% scheme = [scheme1 , scheme2];
% 
% para.name = [scheme,'_ex02_Vesicles_data'];
% pde = ex02_Vesicles_data(para);
% 
% T =20;
% t0 = 0;
% tsave = 0.02*T;
% 
% dt_array = 0.01./2.^(0:4); 
% dt_ref = 0.01/2^3;
% 
% maxIt = length(dt_array);
% 
% % 能量
% lineType = {'b-', 'g--', 'r:'}; 
% 
% figure; % 创建新的图形窗口
% hold on; % 保持当前绘图，以便添加其他元素
% 
% for k = 1:maxIt
%     energy = load([pde.name, '_dt0.00125_S14_energy.txt']);
%     tmp = 1;
% 
%     % 绘制第一条线
%     h1 = plot(energy(tmp:1:end, 1), energy(tmp:1:end, 2), char(lineType(1)), 'LineWidth', 2.5);
%     % 绘制第二条线
%     h2 = plot(energy(tmp:1:end, 1), energy(tmp:1:end, 3), char(lineType(2)), 'LineWidth', 2.5);
%     % 绘制第三条线
%     % h3 = plot(energy(tmp:1:end, 1), energy(tmp:1:end, 4), char(lineType(3)), 'LineWidth', 2.5);
% end
% 
% % 定义图例字符串
% % legend_str = {'Original Energy', 'Modified Energy', 'Discrete modified Energy'};
% % h = legend([h1, h2, h3], legend_str, 'box', 'off', 'Interpreter', 'latex', 'Location', 'north east');
% 
% legend_str = {'Original Energy',  'Discrete Energy'};
% h = legend([h1, h2], legend_str, 'box', 'off', 'Interpreter', 'latex', 'Location', 'north east');
% 
% % 设置标签字体和大小
% xlabel('Time','Fontsize',24);
% ylabel('Energy','Fontsize',24,'interpreter','latex');
% set(gca,'FontSize',12);
% grid on;
% box on;
% set(gca,'linewidth',2.5);
% 
% % 下面部分创建主图
% dirname = '_CAC_Vesicles_data_LM_five_circles';
% dirname = [scheme,dirname];
% datadir = [dirname,'/data'];
% X = load([datadir, '/X.txt']);
% Y = load([datadir, '/Y.txt']);
% 
% %小图1
% t=0;
% ssp = [datadir '/phi_t='  num2str(t) '.txt'];
% phi = load(ssp);
% 
% % 在主图中创建小图
% pos1 = [0.25, 0.6, 0.25, 0.25]; % 位置和大小 [x, y, width, height]
% 
% % 使用 axes 创建小图
% axes('Position', pos1);
% s = pcolor(X, Y, phi);
% s.FaceColor = 'interp';
% s.EdgeColor = 'interp';
% colormap jet;
% axis square;
% axis tight;
% axis off;
% title(['T=',num2str(t)]);
% hold on;
% 
% % 设置其他图形参数
% set(gca, 'FontSize', 12);
% set(gca, 'LineWidth', 1.8);
% 
% annotation('arrow', [0.28, 0.13], [0.73, 0.82]);
% 
% %小图2
% t=2;
% ssp = [datadir '/phi_t='  num2str(t) '.txt'];
% phi = load(ssp);
% 
% % 在主图中创建小图
% pos1 = [0.25, 0.28, 0.25, 0.25]; % 位置和大小 [x, y, width, height]
% 
% % 使用 axes 创建小图
% axes('Position', pos1);
% s = pcolor(X, Y, phi);
% s.FaceColor = 'interp';
% s.EdgeColor = 'interp';
% colormap jet;
% axis square;
% axis tight;
% axis off;
% title(['T=',num2str(t)]);
% hold on;
% 
% % 设置其他图形参数
% set(gca, 'FontSize', 12);
% set(gca, 'LineWidth', 1.8);
% 
% annotation('arrow', [0.28, 0.21], [0.42, 0.3]);
% 
% %小图3
% t=10;
% ssp = [datadir '/phi_t='  num2str(t) '.txt'];
% phi = load(ssp);
% 
% % 在主图中创建小图
% pos1 = [0.47, 0.28, 0.25, 0.25]; % 位置和大小 [x, y, width, height]
% 
% % 使用 axes 创建小图
% axes('Position', pos1);
% s = pcolor(X, Y, phi);
% s.FaceColor = 'interp';
% s.EdgeColor = 'interp';
% colormap jet;
% axis square;
% axis tight;
% axis off;
% title(['T=',num2str(t)]);
% hold on;
% 
% % 设置其他图形参数
% set(gca, 'FontSize', 12);
% set(gca, 'LineWidth', 1.8);
% 
% annotation('arrow', [0.58, 0.52], [0.28, 0.16]);
% 
% %小图4
% t=20;
% ssp = [datadir '/phi_t='  num2str(t) '.txt'];
% phi = load(ssp);
% 
% % 在主图中创建小图
% pos1 = [0.68, 0.28, 0.25, 0.25]; % 位置和大小 [x, y, width, height]
% 
% % 使用 axes 创建小图
% axes('Position', pos1);
% s = pcolor(X, Y, phi);
% s.FaceColor = 'interp';
% s.EdgeColor = 'interp';
% colormap jet;
% axis square;
% axis tight;
% axis off;
% title(['T=',num2str(t)]);
% hold on;
% 
% % 设置其他图形参数
% set(gca, 'FontSize', 12);
% set(gca, 'LineWidth', 1.8);
% 
% annotation('arrow', [0.81, 0.9], [0.28, 0.16]);
% 
% figname = ['D:\paper\phase-field\CAC_Vesicles_surface_LagrangeMultiplier_SAV\fig\',scheme,'N2energy.png'];
% % print(figname,'-dpng', '-r300')
% 
% 
% 
close all;
clear; clc;

% 添加路径
addpath('../','-begin');

% 域参数
domain.left   =   0;
domain.right  =2*pi;
domain.bottom =   0;
domain.top    =2*pi;
Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;

% 网格参数
N = 128;
Nx = N;
Ny = N;

% 物理参数
para.epsilon = 6*pi/128;
para.gamma = 0.1;
para.Re = 1;
para.lambda = 0.1;
para.C0 = 10000;
para.S1 = 4;
para.S2 = 4;
para.S3 = 1;

% 方案设置
scheme1 = 'linear';
% scheme1 = 'nonlinear';
scheme2 = '_bdf2';
scheme = [scheme1, scheme2];
para.name = [scheme,'_ex02_Vesicles_data'];
pde = ex02_Vesicles_data(para);

% 时间参数（仅用于参考，实际未在绘图部分使用）
T = 20;

% 加载能量数据并绘图
energy = load([pde.name, '_dt0.00125_S14_energy.txt']);
lineType = {'b-', 'g--', 'r:'}; 

figure;
hold on;
h1 = plot(energy(:,1), energy(:,2), lineType{1}, 'LineWidth', 2.5);
h2 = plot(energy(:,1), energy(:,3), lineType{2}, 'LineWidth', 2.5);
h3 = plot(energy(:,1), energy(:,4), lineType{3}, 'LineWidth', 2.5);
legend({'Original Energy','Modified Energy','Discrete Energy'}, 'box', 'off', 'Interpreter', 'latex', 'Location', 'north east');
% legend({'Original Energy','Discrete Energy'}, 'box', 'off', 'Interpreter', 'latex', 'Location', 'north east');

xlabel('Time','Fontsize',24);
ylabel('Energy','Fontsize',10,'interpreter','latex');
set(gca,'FontSize',18, 'LineWidth',2.5);


grid on; box on;

% 子图参数设置
timePoints = [0, 2, 10, 20];
positions = [0.2, 0.6, 0.25, 0.25;   % 各子图位置和大小
             0.25, 0.28, 0.25, 0.25;
             0.47, 0.28, 0.25, 0.25;
             0.68, 0.28, 0.25, 0.25];
arrows = {[0.23, 0.13], [0.73, 0.82];  % 箭头坐标
          [0.28, 0.21], [0.42, 0.35];
          [0.58, 0.52], [0.28, 0.2];
          [0.81, 0.90], [0.28, 0.2]};

% 加载几何数据
dirname = [scheme,'_CAC_Vesicles_data_LM_five_circles'];
datadir = fullfile(dirname, 'data');
X = load(fullfile(datadir, 'X.txt'));
Y = load(fullfile(datadir, 'Y.txt'));

% 绘制子图
for i = 1:4
    t = timePoints(i);
    pos = positions(i,:);
    
    % 加载数据
    phiFile = fullfile(datadir, ['phi_t=', num2str(t), '.txt']);
    phi = load(phiFile);
    
    % 创建子图
    axes('Position', pos);
    s = pcolor(X, Y, phi);
    s.FaceColor = 'interp';
    s.EdgeColor = 'none'; % 更清晰的显示
    colormap jet;
    axis square tight off;
    title(['T=', num2str(t)], 'FontSize', 10);
    set(gca, 'LineWidth', 1.8);
    
    % 添加箭头
    annotation('arrow', arrows{i,1}, arrows{i,2}, 'LineWidth', 1.5);
end


figname = ['D:\paper\phase-field\CAC_Vesicles_surface_LagrangeMultiplier_SAV\fig\',scheme,'N2energy.png'];
% print(figname,'-dpng', '-r300')
