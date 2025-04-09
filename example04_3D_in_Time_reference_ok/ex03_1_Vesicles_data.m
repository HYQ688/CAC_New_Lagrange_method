function pde = ex03_1_Vesicles_data(para)
if nargin == 0
    epsilon = 1;
    gamma = 1;
    lambda  = 1;
    C0 = 0;
    Re = 1;
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
    if ~isfield(para,'gamma') || isempty(para.gamma)
        gamma = 1;
    else
        gamma = para.gamma;
    end
    if ~isfield(para,'Re') || isempty(para.Re)
        Re = 1;
    else
        Re = para.Re;
    end 
    if ~isfield(para,'lambda') || isempty(para.lambda)
        lambda = 1;
    else
        lambda = para.lambda;
    end
    if ~isfield(para,'C0') || isempty(para.C0)
        C0 = 0;
    else
        C0 = para.C0;
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
        name = 'ex03_1_Vesicles_data';
    else
        name = para.name;
    end
end

pde = struct('epsilon',epsilon, ...
    'gamma',gamma, ...
    'lambda',lambda, ...
    'Re',Re, ...
    'C0',C0, ...
    'S1',S1, ...
    'S2',S2, ...
    'S3',S3, ...
    'initphi',@initphi, ...
    'initu',@initu,...
    'initv',@initv,...
    'initw',@initw,...
    'initp',@initp,...
    'exactphi',@exactphi, ...
    'exactu',@exactu,...
    'exactv',@exactv,...
    'exactw',@exactw,...
    'exactp',@exactp,...
    'name',name);

    function z = initphi(x,y,z)
        z = exactphi(x,y,z,0);
    end

    function z = exactphi(x,y,z,t)
%         z = (sin(2*x).*sin(2*y)/4+0.48).*(1-sin(t).^2/2);
        z = exp(-t).*sin(2.*x).*cos(2.*y).*cos(2.*z)./4;
    end

    function z = initu(x,y,z)
        z = exactu(x,y,z,0);
    end

    function z = initv(x,y,z)
        z = exactv(x,y,z,0);
    end

    function z = initw(x,y,z)
        z = exactw(x,y,z,0);
    end

    function z = initp(x,y,z)
        z = exactp(x,y,z,0);
    end

    function z = exactu(x,y,z,t)
        z = sin(2*z).*sin(2*y).*sin(x).^2.*sin(t);
    end

    function z = exactv(x,y,z,t)
        z = -sin(2*z).*sin(2*x).*sin(y).^2.*sin(t); 
    end

    function z = exactw(x,y,z,t)
        z = -sin(2*y).*sin(2*x).*sin(z).^2.*sin(t); 
    end

    function z = exactp(x,y,z,t)
        z = cos(x).*sin(t).*sin(y).*sin(z);
    end

   % function z = rhs1(x,y,t,eta2)
   %      z = gamma.*((2.*lambda.*cos(2.*y).*exp(-5.*t).*cos(x).*sin(x).*(exp(4.*t) - 16.*epsilon.^2.*exp(4.*t) + 64.*epsilon.^4.*exp(4.*t) + 48.*cos(2.*y).^4.*cos(x).^4.*sin(x).^4 - 96.*epsilon.^3.*exp(2.*t).*cos(x).^2 + 96.*epsilon.^3.*exp(2.*t).*cos(x).^4 - 24.*epsilon.^3.*cos(2.*y).^2.*exp(2.*t) + 384.*epsilon.^3.*cos(2.*y).^2.*exp(2.*t).*cos(x).^2 - 384.*epsilon.^3.*cos(2.*y).^2.*exp(2.*t).*cos(x).^4 - 16.*cos(2.*y).^2.*exp(2.*t).*cos(x).^2.*sin(x).^2))./epsilon.^3 - (2.*eta2.*lambda.*cos(2.*y).*exp(-3.*t).*cos(x).*sin(x).*(8.*epsilon.^2.*exp(2.*t) - exp(2.*t) + 4.*cos(2.*y).^2.*cos(x).^2.*sin(x).^2))./epsilon) - cos(2.*y).*sin(2.*x).*exp(-t) + 16.*exp(-t).*cos(x).^2.*cos(y).*sin(t).*sin(x).^2.*sin(y).^3 + 4.*exp(-t).*cos(y).*sin(t).*sin(x).^2.*sin(y).*(2.*cos(x).^2 - 1).*(2.*cos(y).^2 - 1);
   % end
   % 
   %  function z = rhs2(x,y,t,eta2)
   %      z = sin(2.*y).*cos(t).*sin(x).^2 - sin(t).*sin(x).*sin(y) + 2.*sin(2.*y).^2.*cos(x).*sin(t).^2.*sin(x).^3 - (2.*sin(2.*y).*cos(x).^2.*sin(t))./Re + (6.*sin(2.*y).*sin(t).*sin(x).^2)./Re - 4.*cos(2.*y).*cos(x).*sin(t).^2.*sin(x).^3.*sin(y).^2 - 96.*lambda.*cos(2.*x).*cos(2.*y).^4.*sin(2.*x).^3.*exp(-4.*t) + 48.*lambda.*cos(2.*x).^3.*cos(2.*y).^4.*sin(2.*x).*exp(-4.*t) + 48.*lambda.*cos(2.*x).*cos(2.*y).^2.*sin(2.*x).^3.*sin(2.*y).^2.*exp(-4.*t) + (8.*lambda.*cos(2.*x).*cos(2.*y).^4.*sin(2.*x).^3.*exp(-4.*t))./epsilon.^3 - (6.*lambda.*cos(2.*x).*cos(2.*y).^6.*sin(2.*x).^5.*exp(-6.*t))./epsilon.^3 - 128.*epsilon.*lambda.*cos(2.*x).*cos(2.*y).^2.*sin(2.*x).*exp(-2.*t) + (32.*lambda.*cos(2.*x).*cos(2.*y).^2.*sin(2.*x).*exp(-2.*t))./epsilon - (2.*lambda.*cos(2.*x).*cos(2.*y).^2.*sin(2.*x).*exp(-2.*t))./epsilon.^3 - (2.*eta2.*lambda.*cos(2.*x).*cos(2.*y).^2.*sin(2.*x).*exp(-2.*t))./epsilon + (2.*eta2.*lambda.*cos(2.*x).*cos(2.*y).^4.*sin(2.*x).^3.*exp(-4.*t))./epsilon + 16.*epsilon.*eta2.*lambda.*cos(2.*x).*cos(2.*y).^2.*sin(2.*x).*exp(-2.*t);
   %  end
   % 
   %  function z = rhs3(x,y,t,eta2)
   %      z = cos(x).*cos(y).*sin(t) - sin(2.*x).*cos(t).*sin(y).^2 + 2.*sin(2.*x).^2.*cos(y).*sin(t).^2.*sin(y).^3 + (2.*sin(2.*x).*cos(y).^2.*sin(t))./Re - (6.*sin(2.*x).*sin(t).*sin(y).^2)./Re - 4.*cos(2.*x).*cos(y).*sin(t).^2.*sin(x).^2.*sin(y).^3 - 48.*lambda.*cos(2.*y).*sin(2.*x).^4.*sin(2.*y).^3.*exp(-4.*t) + 96.*lambda.*cos(2.*y).^3.*sin(2.*x).^4.*sin(2.*y).*exp(-4.*t) - 48.*lambda.*cos(2.*x).^2.*cos(2.*y).^3.*sin(2.*x).^2.*sin(2.*y).*exp(-4.*t) - (8.*lambda.*cos(2.*y).^3.*sin(2.*x).^4.*sin(2.*y).*exp(-4.*t))./epsilon.^3 + (6.*lambda.*cos(2.*y).^5.*sin(2.*x).^6.*sin(2.*y).*exp(-6.*t))./epsilon.^3 + 128.*epsilon.*lambda.*cos(2.*y).*sin(2.*x).^2.*sin(2.*y).*exp(-2.*t) - (32.*lambda.*cos(2.*y).*sin(2.*x).^2.*sin(2.*y).*exp(-2.*t))./epsilon + (2.*lambda.*cos(2.*y).*sin(2.*x).^2.*sin(2.*y).*exp(-2.*t))./epsilon.^3 - 16.*epsilon.*eta2.*lambda.*cos(2.*y).*sin(2.*x).^2.*sin(2.*y).*exp(-2.*t) + (2.*eta2.*lambda.*cos(2.*y).*sin(2.*x).^2.*sin(2.*y).*exp(-2.*t))./epsilon - (2.*eta2.*lambda.*cos(2.*y).^3.*sin(2.*x).^4.*sin(2.*y).*exp(-4.*t))./epsilon;
   %  end



end