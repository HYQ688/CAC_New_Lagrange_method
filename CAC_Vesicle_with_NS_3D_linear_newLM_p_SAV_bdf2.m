function [phi, U, u, v, w, p, Q, eta2] = CAC_Vesicle_with_NS_3D_linear_newLM_p_SAV_bdf2(pde,domain,Nx,Ny,Nz,time,option)

global dt kx ky kz kxx kyy kzz k2 k4 hx hy hz Lx Ly Lz ...
       epsilon gamma Re lambda C0 S1 S2 S3 phi_init

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
    figname_mass = [pde.name, '_dt', num2str(time.dt), '_S1', num2str(pde.S1), '_mass_data.txt'];
    figname_energy = [pde.name, '_dt', num2str(time.dt), '_S1', num2str(pde.S1),'_energy.txt'];
    out1 = fopen(figname_mass,'a');
    out2 = fopen(figname_energy,'a');
end

tol = option.tol;
tolit = option.tolit;
maxit = option.maxit;

%%
T  = time.T;
t  = time.t0;
dt = time.dt;
tsave = time.tsave;

dir_fig  = [pde.name '/fig'];
dir_data = [pde.name '/data'];

epsilon = pde.epsilon;
gamma   = pde.gamma;
lambda  = pde.lambda;
Re      = pde.Re;
C0      = pde.C0;
S1      = pde.S1;
S2      = pde.S2;
S3      = pde.S3;

Lx = domain.xb - domain.xa;
Ly = domain.yb - domain.ya;
Lz = domain.zb - domain.za;

hx = Lx/Nx;
hy = Ly/Ny;
hz = Lz/Nz;

x = linspace( domain.xa, domain.xb-hx, Nx); 
y = linspace( domain.ya, domain.yb-hy, Ny);
z = linspace( domain.za, domain.zb-hz, Nz);

% [k_x,k_y,kx,ky,kxx,kyy,k2,k4] = prepare_fftn_v2(Lx,Ly,Nx,Ny);
k_x = 1i*[0:Nx/2 -Nx/2+1:-1]*(2*pi/Lx);
k_y = 1i*[0:Ny/2 -Ny/2+1:-1]*(2*pi/Ly);
k_z = 1i*[0:Nz/2 -Nz/2+1:-1]*(2*pi/Lz);
[kx, ky, kz] = ndgrid(k_x,k_y,k_z);

k2x = k_x.^2;
k2y = k_y.^2;
k2z = k_z.^2;
[kxx, kyy, kzz] = ndgrid(k2x,k2y,k2z);
k2 = kxx + kyy + kzz;
k4 = k2.^2;

[xx,yy,zz] = ndgrid(x,y,z);
phi0 = pde.initphi(xx,yy,zz);
phi_init = phi0;
u0 = pde.initu(xx,yy,zz);
v0 = pde.initv(xx,yy,zz);
w0 = pde.initw(xx,yy,zz);
p0 = pde.initp(xx,yy,zz);
nfigure =1;

nplot = round((T-t)/dt);
nsave = round(tsave/dt);

tstart = tic;

mu0 = lambda.*( epsilon.*lap_diff(lap_diff(phi0)) + fun_w(phi0));
U0 = U_init(phi0);
Q0  = 1;

%% initialization phi 1Lx.*Ly.*Lz
time1 = time; time1.T = dt;
[phi1,U1,u1,v1,w1,p1,Q1,~] = CAC_Vesicle_with_NS_3D_linear_newLM_p_SAV_1st(pde,domain,Nx,Ny,Nz,time1,option);
mu1 = lambda.*( epsilon.*lap_diff(lap_diff(phi1)) + fun_w(phi1));

t = t+dt;

for nt = 2:nplot
    t = t+dt;
     
    phi_star = 2*phi1-phi0;
    mu_star = 2*mu1-mu0;
    u_star = 2*u1-u0;
    v_star = 2*v1-v0;
    w_star = 2*w1-w0;
    
    
    % step 1
    H = fun_H(phi_star);
    K = (4*U1-U0)/3 - 1/2*fun_inner(H,(4*phi1-phi0)/3);
    
    %step 2 
    if isfield(pde,'rhs1') && isfield(pde,'exactphi')
        rhs1 = pde.rhs1(xx,yy,zz,t,eta_2);
    else
        rhs1 = 0;
    end

    g1 = (4*phi1-phi0)./2./dt./gamma./lambda + rhs1./gamma./lambda  ...
            + S1./epsilon.^3.*(phi_star - 1./(Lx.*Ly.*Lz).*fun_inner(phi_star,1)) ...
            - S2./epsilon.*(lap_diff(phi_star) - 1./(Lx.*Ly.*Lz).*fun_inner(lap_diff(phi_star),1)) ...
            + S3.*epsilon.*(lap_diff(lap_diff(phi_star)) - 1./(Lx.*Ly.*Lz).*fun_inner(lap_diff(lap_diff(phi_star)),1));
    g2 = -(u_star.*diff_x(phi_star) + v_star.*diff_y(phi_star)  + w_star.*diff_z(phi_star))./gamma./lambda;
    g3 = delta_A(phi_star) - 1./(Lx.*Ly.*Lz).*fun_inner(1,delta_A(phi_star));
    g4 = g2;

    psi1 = inv_A(g1-H.*K); psi2 = inv_A(g2);
    psi3 = inv_A(g3);      psi4 = inv_A(g4);

    psih = inv_A(1./2.*H);

    b1 = fun_inner(psi1,H); b2 = fun_inner(psi2,H);
    b3 = fun_inner(psi3,H); b4 = fun_inner(psi4,H);

    a = fun_inner(psih,H);

    X = b1./(1+a); Y = b2./(1+a);
    Z = b3./(1+a); W = b4./(1+a);
    
    %step 3
    phi_11 = psi1 - psih.*X; phi_12 = psi2 - psih.*Y;
    phi_21 = psi3 - psih.*Z; phi_22 = psi4 - psih.*W;

    U_11 = K + 1./2.*X; U_12 = 1./2.*Y;
    U_21 = 1./2.*Z;     U_22 = 1./2.*W;

    mu_11 = get_mu1(phi_11,phi_star,U_11,H);
    mu_12 = get_mu2(phi_12,U_12,H);
    mu_21 = get_mu3(phi_21,phi_star,U_21,H);
    mu_22 = get_mu2(phi_22,U_22,H);
   
   %step 4    
    if isfield(pde,'rhs2') && isfield(pde,'exactu')
        rhs2 = pde.rhs2(xx,yy,zz,t,eta_2);
        rhs3 = pde.rhs3(xx,yy,zz,t,eta_2);
        rhs4 = pde.rhs4(xx,yy,zz,t,eta_2);
    else
        rhs2 = 0;
        rhs3 = 0;
        rhs4 = 0;
    end

    G1_u = rhs2 + (4*u1-u0)./dt/2 - diff_x(p1);  G1_v = rhs3 + (4*v1-v0)./dt/2 - diff_y(p1); G1_w = rhs4 + (4*w1-w0)./dt/2 - diff_z(p1);
    G2_u = -(u_star.*diff_x(u_star) + v_star.*diff_y(u_star) + w_star.*diff_z(u_star)) + mu_star.*diff_x(phi_star);
    G2_v = -(u_star.*diff_x(v_star) + v_star.*diff_y(v_star) + w_star.*diff_z(v_star)) + mu_star.*diff_y(phi_star);
    G2_w = -(u_star.*diff_x(w_star) + v_star.*diff_y(w_star) + w_star.*diff_z(w_star)) + mu_star.*diff_z(phi_star);
    G3_u = 0; G3_v = 0; G3_w = 0;
    G4_u = G2_u; G4_v = G2_v; G4_w = G2_w;

    utilde_11 = inv_B(G1_u); vtilde_11 = inv_B(G1_v); wtilde_11 = inv_B(G1_w);
    utilde_12 = inv_B(G2_u); vtilde_12 = inv_B(G2_v); wtilde_12 = inv_B(G2_w);
    utilde_21 = inv_B(G3_u); vtilde_21 = inv_B(G3_v); wtilde_21 = inv_B(G3_w);
    utilde_22 = inv_B(G4_u); vtilde_22 = inv_B(G4_v); wtilde_22 = inv_B(G4_w);

    %step 5
    nu_11 = get_nu(utilde_11,vtilde_11,wtilde_11,mu_11,phi_star,mu_star,u_star,v_star,w_star) + (4*Q1-Q0)./3;
    nu_12 = get_nu(utilde_12,vtilde_12,wtilde_12,mu_12,phi_star,mu_star,u_star,v_star,w_star);
    Q_1 = nu_11./(1-nu_12);
    
    nu_21 = get_nu(utilde_21,vtilde_21,wtilde_21,mu_21,phi_star,mu_star,u_star,v_star,w_star);
    nu_22 = get_nu(utilde_22,vtilde_22,wtilde_22,mu_22,phi_star,mu_star,u_star,v_star,w_star);
    Q_2 = nu_21./(1-nu_22);

    %stpe 6 
    phi_1 = phi_11 + Q_1*phi_12;  phi_2 = phi_21 + Q_2*phi_22;
    mu_1  = mu_11  + Q_1*mu_12;   mu_2  = mu_21  + Q_2*mu_22;
    U_1   = U_11   + Q_1*U_12;    U_2   = U_21   + Q_2*U_22;
    utilde_1 = utilde_11 + Q_1*utilde_12;vtilde_1 = vtilde_11 + Q_1*vtilde_12;wtilde_1 = wtilde_11 + Q_1*wtilde_12;
    utilde_2 = utilde_21 + Q_2*utilde_22;vtilde_2 = vtilde_21 + Q_2*vtilde_22;wtilde_2 = wtilde_21 + Q_2*wtilde_22;

    %step 7
    eta2 = fun_inner(delta_A(phi_star),4*phi1-phi0-3*phi_1)./ fun_inner(delta_A(phi_star),3*phi_2);
    
    Q  = Q_1 + eta2.*Q_2;     phi = phi_1 + eta2.*phi_2;
    mu = mu_1 + eta2.*mu_2;   U   = U_1 + eta2.*U_2;
    utilde  = utilde_1 + eta2.*utilde_2;  vtilde  = vtilde_1 + eta2.*vtilde_2; wtilde  = wtilde_1 + eta2.*wtilde_2;
    
    
    %step 8
    p1_hat = fftn(p1);
    ut_hat = fftn(utilde);
    vt_hat = fftn(vtilde);
    wt_hat = fftn(wtilde);
    p_hat = p1_hat + (kx.*ut_hat+ky.*vt_hat+kz.*wt_hat)*3/(2*dt)./k2;
    p_hat(1,1,1) = 0;
    p = real(ifftn(p_hat));

    u = utilde -2*dt/3*diff_x(p-p1);
    v = vtilde -2*dt/3*diff_y(p-p1);
    w = wtilde -2*dt/3*diff_z(p-p1);
    
    %% update phi0
    phi0 = phi1; phi1 = phi;
    mu0  = mu1;  mu1  = mu;
    U0   = U1;   U1   = U;
    u0   = u1;   u1   = u;
    v0   = v1;   v1   = v;
    w0   = w1;   w1   = w;
    p0   = p1;   p1   = p;
    Q0   = Q1;   Q1   = Q;
    
    if 1 == option.energyflag
        calculate_energy1(out1,out2,t,u1,u0,v1,v0,w1,w0,p1,phi1,phi0,U1,U0,Q1,Q0);
    end

    if  0 == mod(nt,nsave)
        if 1 == option.printflag
            timeElapsed = toc(tstart);
            fprintf('Q=%.4e,eta2=%.4e,epsilon=%.3f,t=%.5f/%.4f, dt=%.2e, Nx=%d, Ny=%d, timeElapsed=%f\n',Q,eta2,epsilon,t,T,dt,Nx,Ny,timeElapsed);
        end
        
        if 1 == option.saveflag
            if ~exist(dir_data,'dir')
                mkdir(dir_data);
            end
            ss1 = [dir_data '/phi_t=' num2str(t) '.txt'];
            writematrix(phi0,ss1,'Delimiter',' ');
        end
        
        nfigure = nfigure +1;
        if 1 == option.plotflag
            if 1 == option.vtkflag
                write_vtk_grid_values(dir_data,x,y,z,nt,phi0);
            end
            if 1 == option.saveflag
                showsolution_3D(nfigure,xx,yy,zz,phi,t,dir_fig);
            else
                showsolution_3D(nfigure,xx,yy,zz,phi,t);
            end
        end
        if 1 == option.savefinal
        name=[option.scheme,'phi_e',num2str(pde.epsilon),'gamma=',num2str(pde.gamma),...
            'S1=',num2str(pde.S1),'Nx=',num2str(Nx),'Ny=',num2str(Ny),'Nz=',num2str(Nz),'T=',num2str(t)];
        filename=[name '.mat'];
        save(filename,'epsilon','xx','yy','zz','hx','hy','hz','Nx','Ny','Nz','dt','T',...
            'phi','u','v','w','p','domain');
        end
    end
    
end

if 1 == option.savefinal
   name=[option.scheme,'_phi_e',num2str(pde.epsilon),'gamma=',num2str(pde.gamma),'S1=',num2str(pde.S1),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'Nz=',num2str(Nz),'dt=',num2str(dt)];    
    filename=[name '.mat'];
    save(filename,'epsilon','xx','yy','zz','hx','hy','hz','Nx','Ny','Nz','dt','T',...
            'phi','u','v','w','p','domain');
end

if 1 == option.energyflag
    fclose(out1);
    fclose(out2);
end
end


function result = U_init(phi)
global C0
if W(phi) + C0 <0
    disp("Root < 0");
    return;
end
result  = sqrt(W(phi) + C0);
end

function result = fun_H(phi)
global C0 hx hy hz Lx Ly Lz
if W(phi) + C0 <0
    disp("Root < 0");
    return;
end
H  = fun_w(phi)./sqrt(W(phi) + C0);
H_bar = fftn(H);
result = H - 1./(Lx*Ly*Lz).*H_bar(1,1,1)*hx*hy*hz;
end

function result = inv_A(phi)
global dt k2 epsilon gamma lambda Lx Ly Lz S1 S2 S3
    L1 = epsilon.*k2.^2;
    L2 = S1/epsilon.^3 - S2/epsilon.*k2 + S3.*epsilon.*k2.^2;
    phihat = fftn(phi);
    r      = phihat./(3/2/dt/gamma/lambda + L1 + L2);
    r(1,1,1) = phihat(1,1,1)./(3/2/dt/gamma/lambda + L1(1,1,1) + L2(1,1,1) - 1/(Lx*Ly*Lz).*L1(1,1,1).*Lx.*Ly.*Lz - 1/(Lx*Ly*Lz).*L2(1,1,1).*Lx.*Ly.*Lz);
    result = real(ifftn(r));
end

function result = inv_B(u)
global dt Re k2
    uhat = fftn(u);
    r = uhat./(3/2/dt - 1/Re.*k2);
    result = real(ifftn(r));
end


function result = W(phi)
global epsilon
    result = 6./epsilon.^2.*phi.^2.*grad_square(phi) ...
             - 2/epsilon.^2.*grad_square(phi) ...
             + 1/epsilon.^4.*(f(phi)).^2;
    result = result*epsilon/2;
    result = fun_inner(1,result);
end

function result = fun_w(phi)
global epsilon
    div_term = diff_x(phi.^2.*diff_x(phi)) + diff_y(phi.^2.*diff_y(phi)) + diff_z(phi.^2.*diff_z(phi));
    result = 6./epsilon.*(phi.*grad_square(phi) - div_term) ...
             + 2/epsilon.*lap_diff(phi) ...
             + 1/epsilon.^3.*f(phi).*f_der(phi);
end

function result = get_mu1(phi,phi_star,U,H)
    global epsilon lambda S1 S2 S3
    result = lambda.*( epsilon.*lap_diff(lap_diff(phi)) + H.*U + S1./epsilon.^3.*(phi - phi_star) - S2./epsilon.*lap_diff(phi - phi_star) + epsilon.*S3.*lap_diff(lap_diff(phi - phi_star))  );
end

function result = get_mu2(phi,U,H)
    global epsilon lambda S1 S2 S3
    result = lambda.*( epsilon.*lap_diff(lap_diff(phi)) + H.*U + S1./epsilon.^3.*phi - S2./epsilon.*lap_diff(phi) + epsilon.*S3.*lap_diff(lap_diff(phi))  );
end

function result = get_mu3(phi,phi_star,U,H)
    global epsilon lambda S1 S2 S3
    result = lambda.*( epsilon.*lap_diff(lap_diff(phi)) + H.*U - delta_A(phi_star) + S1./epsilon.^3.*phi - S2./epsilon.*lap_diff(phi) + epsilon.*S3.*lap_diff(lap_diff(phi))  );
end

function result = get_nu(u,v,w,mu,phi_star,mu_star,u_star,v_star,w_star)
global dt hx hy hz
   result = fftn((u_star.*diff_x(phi_star)+v_star.*diff_y(phi_star) +w_star.*diff_z(phi_star)).*mu ...
                 -mu_star.*((diff_x(phi_star).*u + diff_y(phi_star).*v + diff_z(phi_star).*w)) ...
                 +(u_star.*diff_x(u_star) + v_star.*diff_y(u_star)  + w_star.*diff_z(u_star)).*u ...
                 +(u_star.*diff_x(v_star) + v_star.*diff_y(v_star)  + w_star.*diff_z(v_star)).*v...
                 +(u_star.*diff_x(w_star) + v_star.*diff_y(w_star)  + w_star.*diff_z(w_star)).*w) ;
    result = dt.*result(1,1,1)*hx*hy*hz;
end

function r = fun_inner(f,g)
global hx hy hz
    r1 = fftn(f.*g);
    r = r1(1,1,1)*hx*hy*hz;
end

function [] = calculate_energy1(out1,out2,t,u1,u0,v1,v0,w1,w0,p1,phi1,phi0,U1,U0,Q1,Q0)
global C0 epsilon lambda S1 S2 S3 dt
rr = 1./2.*epsilon.*fun_inner(1,lap_diff(phi1).^2);
energy_part1 = lambda.*fun_inner(1,1./2.*epsilon.*lap_diff(phi1).^2);

energy_part2 = 1./2.*fun_inner(1,u1.^2 + v1.^2 + w1.^2);

energy_part3 = lambda.*W(phi1);

energy_original = energy_part1 + energy_part2 + energy_part3;

energy_modified = energy_part1 + energy_part2 + lambda*(U1^2-C0) + 1./2*(Q1^2-1);

energy_phi = fun_inner(1,lap_diff(phi1).^2) + fun_inner(1,(2.*lap_diff(phi1)-lap_diff(phi0)).^2);
energy_phi = lambda.*epsilon.*energy_phi./4;

energy_u = fun_inner(1,u1.^2 + v1.^2 + w1.^2) + fun_inner(1,((2*u1-u0).^2+(2*v1-v0).^2+(2*w1-w0).^2));
energy_u = energy_u./4;

energy_p = fun_inner(1,dt.^2/3*grad_square(p1));

energy_U = lambda*(0.5*U1.^2 + 0.5*(2*U1-U0).^2-C0);

energy_Q = 1./2*(0.5*Q1.^2 + 0.5*(2*Q1-Q0).^2) - 0.5;

energy_S = S1./2./epsilon.^3.*(phi1 - phi0).^2 ...
         + S2./2./epsilon.*grad_square(phi1 - phi0) ...
         + S3./2.*epsilon.*lap_diff(phi1 - phi0).^2;
energy_S = lambda*fun_inner(1,energy_S);

energy_discrete = energy_phi + energy_u + energy_p + energy_U + energy_Q + energy_S;
% energy_modified = energy_part1 + energy_u + energy_p + energy_U + energy_Q + energy_S;
% energy_modified = energy_part1 + energy_part2 + lambda*(U1^2 - C0) + 1/2*(Q1^2 - 1);

mass    = fun_inner(1,(phi1 + 1)./2);

surface = A(phi1);

fprintf(out1,'%14.6e  %.8f %.8f \n',t,mass,surface);
fprintf(out2,'%14.6e  %f  %f  %f  %f \n',t,energy_original,energy_modified,energy_discrete,Q1);
end

function lap=lap_diff(phi)
global k2
    lap=real(ifftn((k2.*fftn(phi))));
end

function lap=diff_x(phi)
global kx
    lap=real(ifftn((kx.*fftn(phi))));
end

function lap=diff_y(phi)
global ky
    lap=real(ifftn((ky.*fftn(phi))));
end

function lap=diff_z(phi)
global kz
    lap=real(ifftn((kz.*fftn(phi))));
end

function result = grad_square(phi)
    result = diff_x(phi).^2 + diff_y(phi).^2  + diff_z(phi).^2;
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

function result = A(phi)
global hx hy hz epsilon
    r = fftn(epsilon./2.*grad_square(phi)+1/epsilon.*F(phi));
    result = r(1,1,1)*hx*hy*hz;
end

function result = delta_A(phi)
global epsilon
    result = -epsilon*lap_diff(phi) + 1/epsilon.*f(phi);
end



