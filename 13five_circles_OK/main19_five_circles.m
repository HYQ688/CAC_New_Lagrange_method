close all;
clear; clc;

% add path
addpath('../','-begin');

% Space: Domain and N
% d = 0.54*pi;
domain.left   =   0 ;
domain.right  =2*pi ;
domain.bottom =   0 ;
domain.top    = 2*pi;

Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;

N = 128;
Nx = N;
Ny = N;

% scheme1 = 'linear';
scheme1 = 'nonlinear';

% scheme2 = '_1st';   % First-order scheme
% scheme2 = '_2cn';   % Second-order CN scheme
scheme2 = '_bdf2';  % Second-order BDF scheme

scheme = [scheme1 , scheme2];

% Parameters
para.epsilon = 6*pi/128;
para.gamma = 1;
para.Re = 1;
para.lambda = 0.1;
para.C0 = 10000;
para.S1 = 4;
para.S2 = 4;
para.S3 = 1;

para.name = [scheme,'_CAC_Vesicles_data_LM_six_circles'];


% Time: dt T
% T = 0.2;
% t0 = 0;
% tsave = 0.02*T;
% 
% dt_array = 1./2.^(0:4); 
% dt_ref = 1e-6;

T =20;
t0 = 0;
tsave = 0.02*T;

dt_array = 0.01./2.^(0:4); 
dt_ref = 0.01/2^3;

maxIt = length(dt_array);

%% option
option.scheme = scheme;
option.plotflag  = 0;
option.printflag = 1;
option.vtkflag  = 0;
option.saveflag  = 1;
option.savefinal  = 1;
option.energyflag = 1;
option.tol = 1e-14;
option.tolit = 1e-11;
option.maxit = 2000;

pde = ex02_Vesicles_data(para);

if 1 == strcmp(scheme,'linear_1st')
    solver_fun = @CAC_Vesicle_2D_linear_newLM_p_SAV_1st;
elseif 1 == strcmp(scheme,'linear_2cn')
    solver_fun = @CAC_Vesicle_2D_linear_newLM_p_SAV_2cn;
elseif 1 == strcmp(scheme,'linear_bdf2')
    solver_fun = @CAC_Vesicle_with_NS_2D_linear_newLM_p_SAV_bdf2;
elseif 1 == strcmp(scheme,'nonlinear_1st')
    solver_fun = @CAC_Vesicle_2D_newLM_p_SAV_1st;
elseif 1 == strcmp(scheme,'nonlinear_2cn')
    solver_fun = @CAC_Vesicle_2D_newLM_p_SAV_2cn;
elseif 1 == strcmp(scheme,'nonlinear_bdf2')
    solver_fun = @CAC_Vesicle_with_NS_2D_newLM_p_SAV_bdf2;
end


%% Run:
if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
    time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
    solver_fun(pde,domain,Nx,Ny,time,option);
end


