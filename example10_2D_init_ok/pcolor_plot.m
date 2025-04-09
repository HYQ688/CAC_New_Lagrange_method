clear;
clc;

scheme1 = 'linear';
% scheme1 = 'nonlinear';

% scheme2 = '_1st';   % First-order scheme
% scheme2 = '_2cn';   % Second-order CN scheme
scheme2 = '_bdf2';  % Second-order BDF scheme

scheme = [scheme1 , scheme2];

dirname = '_CAC_Vesicles_data_LM';
dirname = [scheme,dirname];

datadir = [dirname,'/data'];
figdir  = [dirname,'/',dirname];

X = load([datadir, '/X.txt']);
Y = load([datadir, '/Y.txt']);

figure(6)
for t= [0 1.2 4 10 20]
    t
    ssp = [datadir '/phi_t=' num2str(t) '.txt'];
    phi = load(ssp);

%     subplot(1,2,1)
    s=pcolor(X,Y,phi);
    s.FaceColor='interp';
    s.EdgeColor='interp';
%     view(0,90);
    colormap jet;
%     colormap viridis;
    axis square;
    axis tight;
    axis off;
%     subplot(1,2,2)
%     s=pcolor(lap_diff(phi,k2));
%     s.FaceColor='interp';
%     s.EdgeColor='interp';
% %     view(0,90);
%     colormap gray;
%     axis square;
%     axis tight;
%     axis off;
    % colorbar('Position',[0.845 0.18 0.03 0.66],'Fontsize',15);
    % caxis([-1 1])
%     colorbar off;
    
    drawnow;
    pause(0.25)
    
%     figname = [figdir '/phi_t=' num2str(t) '.eps'];
%     print(figname,'-depsc2', '-r120')

    figname = ['D:\paper\phase-field\CAC_Vesicles_surface_LagrangeMultiplier_SAV\fig\',dirname '_phi_t=' num2str(t) '.png'];
    % print(figname,'-dpng', '-r300')
end


