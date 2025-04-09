function pde = ex03_Vesicles_data(para)
if nargin == 0
    epsilon = 1;
    M = 1;
    C0 = 0;
    beta_m = 1;
    e1  = 0;
    e2 = 0;
    S1 = 0;
    S2 = 0;
    S3 = 0;
    name = '';
else
    if ~isstruct(para)
        exit('we need a struct data');
    end
    if ~isfield(para,'epsilon') || isempty(para.epsilon)
        epsilon = 1;
    else
        epsilon = para.epsilon;
    end
    if ~isfield(para,'M') || isempty(para.M)
        M = 1;
    else
        M = para.M;
    end
    if ~isfield(para,'gamma') || isempty(para.gamma)
        gamma = 1;
    else
        gamma = para.gamma;
    end
    if ~isfield(para,'lambda') || isempty(para.lambda)
        lambda = 1;
    else
        lambda = para.lambda;
    end
    if ~isfield(para,'nu') || isempty(para.nu)
        nu = 0;
    else
        nu = para.nu;
    end
    if ~isfield(para,'C0') || isempty(para.C0)
        C0 = 0;
    else
        C0 = para.C0;
    end
    if ~isfield(para,'beta_m') || isempty(para.beta_m)
        beta_m = 1;
    else
        beta_m = para.beta_m;
    end
    if ~isfield(para,'e1') || isempty(para.e1)
        e1 = 0;
    else
        e1 = para.e1;
    end
    if ~isfield(para,'e2') || isempty(para.e2)
        e2 = 0;
    else
        e2 = para.e2;
    end
    if ~isfield(para,'S1') || isempty(para.S1)
        S1 = 0;
    else
        S1 = para.S1;
    end
    if ~isfield(para,'S2') || isempty(para.S2)
        S2 = 0;
    else
        S2 = para.S2;
    end
    if ~isfield(para,'S3') || isempty(para.S3)
        S3 = 0;
    else
        S3 = para.S3;
    end
    if ~isfield(para,'name') || isempty(para.name)
        name = 'ex02_Vesicles_data';
    else
        name = para.name;
    end
end

pde = struct('epsilon',epsilon, ...
    'M',M, ...
    'C0',C0, ...
    'nu',nu,...
    'beta_m', beta_m, ... 
    'gamma',gamma,...
    'lambda',lambda,...
    'e1',e1, ...
    'e2',e2, ...
    'S1',S1, ...
    'S2',S2, ...
    'S3',S3, ...
    'initphi',@initphi, ...
    'initu',@initu,...
    'initv',@initv,...
    'initp',@initp,...
    'name',name);

   function z = initphi(x,y)
        d = 0.18*pi;
        rl =[0.22*pi,0.22*pi,0.22*pi,0.22*pi,0.22*pi,0.22*pi];
        xl = [pi,pi+3*d,pi+3*cos(2*pi/5)*d,pi+3*cos(4*pi/5)*d,pi+3*cos(6*pi/5)*d,pi+3*cos(8*pi/5)*d];
        yl = [pi,pi,pi+3*sin(2*pi/5)*d,pi+3*sin(4*pi/5)*d,pi+3*sin(6*pi/5)*d,pi+3*sin(8*pi/5)*d];
        rx = 1;
        ry = 1;
        
        z = tanh((rl(2)-sqrt((x-xl(2)).^2./rx.^2+(y-yl(2)).^2./ry.^2))./(sqrt(2)*epsilon)) ...
            +tanh((rl(3)-sqrt((x-xl(3)).^2./rx.^2+(y-yl(3)).^2./ry.^2))./(sqrt(2)*epsilon)) ...
            +tanh((rl(4)-sqrt((x-xl(4)).^2./rx.^2+(y-yl(4)).^2./ry.^2))./(sqrt(2)*epsilon)) ...
            +tanh((rl(5)-sqrt((x-xl(5)).^2./rx.^2+(y-yl(5)).^2./ry.^2))./(sqrt(2)*epsilon)) ...
            +tanh((rl(6)-sqrt((x-xl(6)).^2./rx.^2+(y-yl(6)).^2./ry.^2))./(sqrt(2)*epsilon))+4;
    end
    
    function z = initu(x,y)
        z = zeros(size(x));
    end

    function z = initv(x,y)
        z = zeros(size(x));
    end

    function z = initp(x,y)
        z = zeros(size(x));
    end

end