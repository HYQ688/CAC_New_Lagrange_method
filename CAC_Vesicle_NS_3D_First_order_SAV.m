function [phi, mu, U, u, v, p, Q] = CAC_Vesicle_NS_3D_First_order_SAV(pde,domain,Nx,Ny,Nz,time,option)
% Solve 2D Vesicle model
% Qi Li
% Reference: Yang, Xiaofeng. 
% "Numerical approximations of the Navier–Stokes equation coupled with volume-conserved multi-phase-field vesicles system: fully-decoupled, linear, unconditionally energy stable and second-order time-accurate numerical scheme." 
% Computer Methods in Applied Mechanics and Engineering 375 (2021): 113600.
% 04/30/2023
global dt kx ky kz kxx kyy kzz k2 k4 hx hy hz Lx Ly Lz epsilon ... % gamma  
       beta ...
       C0 nu lambda M e1 e2 S1 S2 S3

if ~exist('option','var'), option = []; end
if ~isfield(option,'tol')
    option.tol = 10^-14;   % default tol
end
if ~isfield(option,'tolit')
    option.tolit = 10^-8;   % default tolit
end
if ~isfield(option,'maxit')
    option.maxit = 2000;   % default maxit
end
if ~isfield(option,'plotflag')
    option.plotflag = 0;   
end
if ~isfield(option,'saveflag')
    option.saveflag = 0;  
end
if ~isfield(option,'savefinal')
    option.savefinal = 0;  
end
if ~isfield(option,'printflag')
    option.printflag = 0;   
end
if ~isfield(option,'vtkflag')
    option.printflag = 0;   
end
if ~isfield(option,'energyflag')
    option.energyflag = 0;   
end
if 1 == option.energyflag
    figname_mass = [pde.name,'_dt_',num2str(time.dt),'_mass.txt'];
    figname_energy = [pde.name,'_dt_',num2str(time.dt),'_energy.txt'];      
    out1 = fopen(figname_mass,'w');
    out2 = fopen(figname_energy,'w');
end

% tol = option.tol;
% tolit = option.tolit;
% maxit = option.maxit;

%%
T  = time.T;
t  = time.t0;
dt = time.dt;
tsave = time.tsave;

dir_fig  = [pde.name '/fig'];
dir_data = [pde.name '/data'];

epsilon = pde.epsilon;
gamma   = pde.gamma;
C0      = pde.C0;
beta_m  = pde.beta_m;
lambda  = pde.lambda;
nu   = pde.nu;
M    = pde.M;
e1   = pde.e1;
e2   = pde.e2;
S1   = pde.S1;
S2   = pde.S2;
S3   = pde.S3;

Lx = domain.xb - domain.xa;
Ly = domain.yb - domain.ya;
Lz = domain.zb - domain.za;

hx = Lx/Nx;
hy = Ly/Ny;
hz = Lz/Nz;

x = linspace( domain.xa, domain.xb-hx, Nx); 
y = linspace( domain.ya, domain.yb-hy, Ny);
z = linspace( domain.za, domain.zb-hz, Nz);

% [k_x,k_y,kx,ky,kxx,kyy,k2,k4] = prepare_fft2_v2(Lx,Ly,Nx,Ny);

k_x = 1i*[0:Nx/2 -Nx/2+1:-1]*(2*pi/Lx);
k_y = 1i*[0:Ny/2 -Ny/2+1:-1]*(2*pi/Ly);
k_z = 1i*[0:Nz/2 -Nz/2+1:-1]*(2*pi/Lz);
[kx, ky] = ndgrid(k_x,k_y);

k2x = k_x.^2;
k2y = k_y.^2;
[kxx, kyy] = ndgrid(k2x,k2y);
k2 = kxx + kyy;
k4 = k2.^2;

[xx,yy] = ndgrid(x,y);
phi0 = pde.initphi(xx,yy);
phi_init = phi0;
u0 = pde.initu(xx,yy);
v0 = pde.initv(xx,yy);
p0 = pde.initp(xx,yy);
nfigure =1;

beta = beta_m*Area(phi0);

%% plot initial value
if 1 == option.saveflag
    if ~exist(dir_data,'dir')
        mkdir(dir_data);
    end
    ss = [dir_data '/phi_t=' num2str(t) '.txt'];
    fid = fopen(ss, 'wt');
    fprintf(fid, '%f\n', phi0(:));
    fclose(fid);
end
if 1 == option.plotflag
    if 1 == option.saveflag
        showsolution_2D(nfigure,xx,yy,phi0,t,dir_fig);
    else
        showsolution_2D(nfigure,xx,yy,phi0,t);
    end
end

nplot = round((T-t)/dt);
nsave = round(tsave/dt);

tstart = tic;

mu0 = -epsilon.*lap_diff(omega(phi0)) + 1./epsilon.*f_der(phi0).*omega(phi0) ...
     + epsilon.*M.*(Area(phi0)-beta).*omega(phi0); 
U0  = fun_U_init(phi0);
Q0  = 1;

% Initial energy
if 1 == option.energyflag
    calculate_energy1(out1,out2,hx,hy,t,u0,v0,p0,phi0,U0,Q0);
end

for nt = 1:nplot
    t = t+dt;
    phi_star = phi0;
    mu_star  = mu0;
    u_star   = u0;
    v_star   = v0;
    
    % Stage 1
    H = fun_H(phi_star);
    G1 = get_rhs1(phi0,gamma);
    G2 = get_rhs2(H);
    G3 = get_rhs3(u0,v0,phi0,gamma);
%     G4 = G2;
   if isfield(pde,'rhs1') && isfield(pde,'exactphi')
        tmp1 = pde.rhs1(xx,yy,t,beta);
        tmp2 = pde.rhs2(xx,yy,t,beta);
        tmp3 = pde.rhs3(xx,yy,t,beta);
        G1 = G1 + tmp1/gamma;
%         rhs2 = rhs2 + tmp2/gamma;
%         rhs3 = rhs3 + tmp3/gamma;
   else
        tmp2 = 0;
        tmp3 = 0;
   end

    phi11 = inv_A(G1,gamma);
    phi12 = inv_A(G2,gamma);
    phi21 = inv_A(G3,gamma);
    phi22 = phi12; % phi22 = inv_A(G4); G2=G4;

    mu11 = epsilon * e1 * lap_diff(lap_diff(phi11)) ...
           + epsilon * e2 * phi11 ...
           + S1/epsilon.^3.*(phi11-phi_star) ...
           - S2/epsilon*lap_diff(phi11-phi_star) ...
           + S3*epsilon*lap_diff(lap_diff(phi11-phi_star));
    mu12 =  epsilon * e1 * lap_diff(lap_diff(phi12)) ...
           + epsilon * e2 * phi12 ...
           + S1/epsilon.^3.*phi12 ...
           - S2/epsilon*lap_diff(phi12) ...
           + S3*epsilon*lap_diff(lap_diff(phi12)) ...
           + epsilon*H;
    mu21 = epsilon * e1 * lap_diff(lap_diff(phi21)) ...
           + epsilon * e2 * phi21 ...
           + S1/epsilon.^3.*phi21 ...
           - S2/epsilon*lap_diff(phi21) ...
           + S3*epsilon*lap_diff(lap_diff(phi21));
    mu22 = epsilon * e1 * lap_diff(lap_diff(phi22)) ...
           + epsilon * e2 * phi22 ...
           + S1/epsilon.^3.*phi22 ...
           - S2/epsilon*lap_diff(phi22) ...
           + S3*epsilon*lap_diff(lap_diff(phi22)) ...
           + epsilon*H;

    % Stage 2
    U_new1 = fun_U1(phi0,phi11,phi12,U0,H);
    U_new2 = fun_U2(phi21,phi22,H);

    % Stage 3
    phi_new1 = phi11 + U_new1*phi12;
    phi_new2 = phi21 + U_new2*phi22;
    mu_new1  = mu11  + U_new1*mu12;
    mu_new2  = mu21  + U_new2*mu22;

    % Stage 4
    [utilde_new1,vtilde_new1] = get_utilde1(p0,u0,v0,tmp2,tmp3);
    [utilde_new2,vtilde_new2] = get_utilde2(u_star,v_star,phi_star,mu_star);

    % Stage 5
    nu1 = fft2(lambda.*(u_star.*diff_x(phi_star)+v_star.*diff_y(phi_star)).*mu_new1...
              -lambda.*mu_star.*((diff_x(phi_star).*utilde_new1 + diff_y(phi_star).*vtilde_new1))...
              +(u_star.*diff_x(u_star) + v_star.*diff_y(u_star)).*utilde_new1 ...
              +(u_star.*diff_x(v_star) + v_star.*diff_y(v_star)).*vtilde_new1);
    nu1 = nu1(1,1)*hx*hy;
    nu2 = fft2(lambda.*(u_star.*diff_x(phi_star)+v_star.*diff_y(phi_star)).*mu_new2...
              -lambda.*mu_star.*((diff_x(phi_star).*utilde_new2 + diff_y(phi_star).*vtilde_new2))...
              +(u_star.*diff_x(u_star) + v_star.*diff_y(u_star)).*utilde_new2 ...
              +(u_star.*diff_x(v_star) + v_star.*diff_y(v_star)).*vtilde_new2);
    nu2 = nu2(1,1)*hx*hy;
    Q = (Q0./dt+nu1)./(1/dt-nu2);

    % Stage 6
    phi = phi_new1 + Q.*phi_new2;
    mu  = mu_new1 + Q.*mu_new2;
    U   = U_new1+ Q.*U_new2;
    utilde = utilde_new1 + Q.*utilde_new2;
    vtilde = vtilde_new1 + Q.*vtilde_new2;
    
    % Stage 7
    p0_hat = fft2(p0);
    ut_hat = fft2(utilde);
    vt_hat = fft2(vtilde);
    p_hat = p0_hat + (kx.*ut_hat+ky.*vt_hat)./k2/dt;
    p_hat(1,1) = 0;
    p = real(ifft2(p_hat));

    u = utilde -dt*diff_x(p-p0);
    v = vtilde -dt*diff_y(p-p0);

    %% update phi0
    phi0 = phi;
    mu0  = mu;
    U0   = U;
    u0   = u;
    v0   = v;
    p0   = p;
    Q0   = Q;

    if 1 == option.energyflag
        calculate_energy1(out1,out2,hx,hy,t,u,v,p,phi,U,Q);
    end

    if  0 == mod(nt,nsave)
        if 1 == option.printflag
            timeElapsed = toc(tstart);
            fprintf('epsilon=%.3f,t=%.4f/%.3f, dt=%.2e, Nx=%d, Ny=%d, timeElapsed=%f\n',epsilon,t,T,dt,Nx,Ny,timeElapsed);
        end
        
        if 1 == option.saveflag
            ss = [dir_data '/phi_t=' num2str(t) '.txt'];
            fid = fopen(ss, 'wt');
            fprintf(fid, '%f\n', phi(:));
            fclose(fid);
        end
        
        nfigure = nfigure +1;
        if 1 == option.plotflag
            if 1 == option.vtkflag
                write_vtk_grid_values(dir_data,x,y,nt,phi0);
            end
            if 1 == option.saveflag
                showsolution_2D(nfigure,xx,yy,phi,t,dir_fig);
            else
%                 subplot(1,2,1);
%                 phiphi = pde.exactphi(xx,yy,t);
                showsolution_2D(nfigure,xx,yy,phi,t);
%                 subplot(1,2,2);
%                 showsolution_2D(nfigure,xx,yy,pde.exact(xx,yy,t),t);
            end
        end
    end
    
end

if 1 == option.savefinal
    name=['phi_e',num2str(pde.epsilon),'gamma',num2str(pde.gamma),...
          'M',num2str(pde.M),'Nx=',num2str(Nx),'Ny=',num2str(Ny),...
          'dt=',num2str(dt),'T=',num2str(T)];    
    filename=[name '.mat'];
    save(filename,'epsilon','xx','yy','hx','hy','Nx','Ny','dt','T',...
         'phi','u','v','p','domain');
end

if 1 == option.energyflag
    fclose(out1);
    fclose(out2);
end
end

function result = fun_U_init(phi)
global hx hy epsilon beta M C0
E1 = fft2(W_tilde(phi));
if E1(1,1)*hx*hy + M./(2*epsilon)*(Area(phi) - beta).^2 + C0 <0
    disp("Root < 0");
    return;
end
result  = sqrt(E1(1,1)*hx*hy + M./(2*epsilon)*(Area(phi) - beta).^2 + C0);
end

function result = W(phi)
global epsilon 
    result = 1./2*(lap_diff(phi)-1./epsilon.^2*f(phi)).^2;
end

function result = W_tilde(phi)
global e1 e2
    result = W(phi) - (e1./2*lap_diff(phi).^2 + e2./2*phi.^2);
end

function result = fun_H(phi)
global epsilon beta C0 M e1 e2 hx hy
E1 = fft2(W_tilde(phi));
if E1(1,1)*hx*hy + M./(2*epsilon)*(Area(phi) - beta).^2 + C0 <0
    disp("Root < 0");
    return;
end
omega_tilde = -e1*lap_diff(lap_diff(phi)) - e2*phi ...
            -lap_diff(omega(phi)) + 1/epsilon.^2*f_der(phi).*omega(phi) ...
            + M*(Area(phi) - beta).*(-lap_diff(phi) + (1/epsilon.^2).*f(phi));
result = omega_tilde./sqrt(E1(1,1)*hx*hy + M./(2*epsilon)*(Area(phi) - beta).^2 + C0);
end

function result = fun_U1(phi0,phi11,phi12,U0,H)
global hx hy
Hphi0  = fft2(H.*phi0);
Hphi0  = Hphi0(1,1)*hx*hy;
Hphi11 = fft2(H.*phi11);
Hphi11 = Hphi11(1,1)*hx*hy;
Hphi12 = fft2(H.*phi12);
Hphi12 = Hphi12(1,1)*hx*hy;
G5 = U0 - 1./2*Hphi0;
result  = (1./2*Hphi11 + G5)./(1 - 1./2*Hphi12);
end

function result = fun_U2(phi21,phi22,H)
global hx hy
Hphi21 = fft2(H.*phi21);
Hphi21 = Hphi21(1,1)*hx*hy;
Hphi22 = fft2(H.*phi22);
Hphi22 = Hphi22(1,1)*hx*hy;
result = 1./2*Hphi21./(1 - 1./2*Hphi22);
end

function [u,v] = get_utilde1(p0,u0,v0,tmp2,tmp3)
global dt k2 nu
rhs1 = -diff_x(p0) + u0/dt+tmp2;
rhs2 = -diff_y(p0) + v0/dt+tmp3;
u = fft2(rhs1)./(1/dt - nu*k2);
v = fft2(rhs2)./(1/dt - nu*k2);
u = real(ifft2(u));
v = real(ifft2(v));
end

function [u,v] = get_utilde2(u_star,v_star,phi_star,mu_star)
global dt k2 nu lambda
rhs1 = -(u_star.*diff_x(u_star) + v_star.*diff_y(u_star)) + lambda*mu_star.*diff_x(phi_star);
rhs2 = -(u_star.*diff_x(v_star) + v_star.*diff_y(v_star)) + lambda*mu_star.*diff_y(phi_star);
u = fft2(rhs1)./(1/dt -nu*k2);
v = fft2(rhs2)./(1/dt -nu*k2);
u = real(ifft2(u));
v = real(ifft2(v));
end

function result = get_rhs1(phi0,gamma)
global dt epsilon e2 S1 S2 S3 hx hy Lx Ly
phi0bar = fft2(phi0);
result = phi0/dt/gamma ...
         + S1/epsilon.^3.*phi0 - S2/epsilon*lap_diff(phi0) ...
         + S3*epsilon*lap_diff(lap_diff(phi0)) ...
         + epsilon*e2*phi0bar(1,1)*hx*hy/(Lx*Ly);
end

function result = get_rhs2(H)
global epsilon hx hy Lx Ly
Hbar = fft2(H);
Hb  = H - Hbar(1,1)*hx*hy/(Lx*Ly);
result = -epsilon*Hb;
end

function result = get_rhs3(u_star,v_star,phi_star,gamma)
result = -(u_star.*diff_x(phi_star) + v_star.*diff_y(phi_star));
result = result/gamma;
end

function result = inv_A(phi,gamma)
global dt k2 epsilon e1 e2 S1 S2 S3
    L1 = epsilon * e1 * k2.^2 + epsilon * e2;
    L2 = S1/epsilon.^3 - S2/epsilon*k2 + S3*epsilon*k2.^2;
    phihat = fft2(phi);
    r      = phihat./(1/dt/gamma + L1 + L2);
    result = real(ifft2(r));
end

function [] = calculate_energy1(out1,out2,hx,hy,t,u1,v1,p1,phi1,U1,Q)
global dt epsilon lambda e1 e2 C0
energye1e2 = fft2(e1/2*lap_diff(phi1).^2 + e2/2*phi1.^2);
energye1e2 = energye1e2(1,1)*hx*hy;
energy_u_1 = 1/2*fft2(u1.^2 + v1.^2);
energy_u = energy_u_1(1,1)*hx*hy;
energy_p = fft2(dt.^2/2*grad_square(p1));
energy_p = energy_p(1,1)*hx*hy;
energy_U = lambda*epsilon*(U1.^2-C0);
energy_Q = 1./2*Q.^2;

energy_modified = energye1e2 + energy_u + energy_p + energy_U + energy_Q;

energy_u = 1/2*fft2(u1.^2 + v1.^2);
energy_W = fft2(W(phi1));
energy_original = energy_u(1,1)*hx*hy + lambda*epsilon*(energy_W(1,1)*hx*hy);

mass = fft2(phi1);
mass = mass(1,1)*hx*hy;

surface = Area(phi1);

fprintf(out1,'%14.6e  %.8f  %.8f\n',t,mass,surface);
fprintf(out2,'%14.6e  %f  %f\n',t,energy_original,energy_modified);
end

function lap=lap_diff(phi)
global k2
    lap=real(ifft2((k2.*fft2(phi))));
end

function lap=diff_x(phi)
global kx
    lap=real(ifft2((kx.*fft2(phi))));
end

function lap=diff_y(phi)
global ky
    lap=real(ifft2((ky.*fft2(phi))));
end

function result = grad_square(phi)
    result = diff_x(phi).^2 + diff_y(phi).^2;
end

function result = omega(phi)
global epsilon
    result = -lap_diff(phi) + 1/(epsilon.^2).*f(phi);
end

function result = f_der(phi)
    result = 3.*phi.^2-1;
end

function result = f(phi)
    result = phi.^3 - phi;
end

function result = F(phi)
    result = 1/4*(phi.^2-1).^2;
end

function result = Area(phi)
global hx hy epsilon
    r = fft2(epsilon./2.*grad_square(phi)+1/epsilon.*F(phi));
    result = r(1,1)*hx*hy;
end
