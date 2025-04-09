clear;
close all; clc;


filename1 = 'nonlinear_bdf2phi_e0.14726gamma=1S1=4Nx=64Ny=64';



for time = 0:0.4:20

filename = [filename1,'T=',num2str(time)];

Nx = 64;
Ny = Nx;
domain.left   =  0;
domain.right  = 2*pi;
domain.bottom =  0;
domain.top    = 2*pi;
Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;
hx = Lx/Nx;
hy = Ly/Ny;

x  = domain.left   + hx*(0:Nx-1);
y  = domain.bottom + hy*(0:Ny-1);

[xx,yy] = ndgrid(x,y);

load(filename,'phi','u','v','p','hx','hy');

figure(1)
contour(xx,yy,phi,[0,0],'r','linewidth',2);


temp = 4;
quiver(xx(1:temp:end,1:temp:end),yy(1:temp:end,1:temp:end),...
       u(1:temp:end,1:temp:end),v(1:temp:end,1:temp:end),...
       1.8,'linewidth',1.2)

% axis([0.65 5.65 0.65 5.65]);
axis square;
box on;
% axis off;
set(gca,'xtick',[],'ytick',[]);

figname1 = ['C:\Users\Administrator\Desktop\论文\figs\速度场\other\',filename,'_uv_a','.png'];
% print(figname1,'-dpng', '-r300')


figure(2)
s=pcolor(xx,yy,p);
s.FaceColor='interp';
s.EdgeColor='interp';
colormap jet;
% colorbar;
colorbar('Position',[0.847 0.16 0.03 0.70],'Fontsize',20);
% colormap viridis;
% colormap parula;
axis square;
box on;
% axis tight;
% axis off;
set(gca,'xtick',[],'ytick',[]);

figname2 = ['C:\Users\Administrator\Desktop\论文\figs\压力\other\',filename,'_p_a','.png'];
% print(figname2,'-dpng', '-r300')

end