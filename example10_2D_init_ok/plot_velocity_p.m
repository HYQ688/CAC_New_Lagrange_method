clear;
close all; clc;


filename1 = 'nonlinear_bdf2phi_e0.14726gamma=1S1=4Nx=128Ny=128';



for time = [0.4]

filename = [filename1,'T=',num2str(time)];
Nx = 128;
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
contour(xx,yy,phi,[-3*pi/128,3*pi/128],'r','linewidth',2);
hold on


temp = 4;
quiver(xx(1:temp:end,1:temp:end),yy(1:temp:end,1:temp:end),...
       u(1:temp:end,1:temp:end),v(1:temp:end,1:temp:end),...
       1.8,'linewidth',1.2)

% axis([0.65 5.65 0.65 5.65]);
axis square;
box on;
% axis off;
set(gca,'xtick',[],'ytick',[]);

figname1 = ['D:\paper\phase-field\CAC_Vesicles_surface_LagrangeMultiplier_SAV\fig\',filename,'_uv_a','.png'];
% print(figname1,'-dpng', '-r300')


figure(3)
s=pcolor(xx,yy,p);
s.FaceColor='interp';
s.EdgeColor='interp';
colormap jet;
% colorbar;
% caxis([-1e-2, 2e-2]);
colorbar('Position',[0.847 0.16 0.03 0.70],'Fontsize',20);

% colormap viridis;
% colormap parula;
axis square;
box on;
% axis tight;
% axis off;
set(gca,'xtick',[],'ytick',[]);

figname2 = ['D:\paper\phase-field\CAC_Vesicles_surface_LagrangeMultiplier_SAV\fig\',filename,'_p_a','.png'];
% print(figname2,'-dpng', '-r300')
% 
end