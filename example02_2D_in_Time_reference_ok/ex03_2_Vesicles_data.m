function pde = ex03_2_Vesicles_data(para)
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
    'initp',@initp,...
    'exactphi',@exactphi, ...
    'exactu',@exactu,...
    'exactv',@exactv,...
    'exactp',@exactp,...
    'rhs_u',@rhs_u,...
    'rhs_v',@rhs_v,...
    'name',name);


    function z = initphi(x,y)
        rl =[0.28*pi,0.28*pi];
        xl = [pi,pi];
        yl = [1.29*pi, -0.29*pi+pi];
        
        z =  tanh((rl(1)-sqrt( (x-xl(1)).^2+(y-yl(1)).^2 ))./(sqrt(2)*epsilon)) ...
            +tanh((rl(2)-sqrt( (x-xl(2)).^2+(y-yl(2)).^2 ))./(sqrt(2)*epsilon))+1;
    end

    function z = initu(x,y)
        z = exactu(x,y);
    end

    function z = initv(x,y)
        z = exactv(x,y);
    end

    function z = initp(x,y)
        z = exactp(x,y);
    end

    % function z = exactu(x,y,t)
    %     z = sin(2*y).*sin(x).^2.*sin(t);
    % end
    % 
    % function z = exactv(x,y,t)
    %     z = -sin(2*x).*sin(y).^2.*sin(t); 
    % end
    % 
    % function z = exactp(x,y,t)
    %     z = cos(x).*sin(t).*sin(y);
    % end

    function z = exactu(x,y)
        z = zeros(size(x));
    end

    function z = exactv(x,y)
        z = zeros(size(x));
    end

    function z = exactp(x,y)
        z = zeros(size(x));
    end

end



