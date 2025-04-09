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

N = 64;
Nx = N;
Ny = N;
Nz = N;

% Parameters
para.epsilon = 6*pi/128;
% para.epsilon = 0.08;
para.gamma = 1;
para.Re = 1;
% para.lambda = 0.001;
para.lambda = 0.1;
para.S1 = 4;
para.S2 = 4;
para.S3 = 1;
para.C0 = 10000;
para.name = 'linear_ex20_BDF2_3D_CAC_Vesicles_data_LM';

% Time: dt T
T = 20;
t0 = 0;
tsave = 0.02*T;

dt_array = 0.01./2.^(0:1); 
dt_ref = 0.01/2^3;

maxIt = length(dt_array);

%% option
option.plotflag  = 0;
option.printflag = 1;
option.vtkflag  = 0;
option.saveflag  = 1;
option.savefinal  = 0;
option.energyflag = 0;
option.tol = 1e-14;
option.tolit = 1e-11;
option.maxit = 2000;

pde = ex02_Vesicles_data(para);

%% Run:
if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
    time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
    CAC_Vesicle_with_NS_3D_linear_newLM_p_SAV_bdf2(pde,domain,Nx,Ny,Nz,time,option);
end

