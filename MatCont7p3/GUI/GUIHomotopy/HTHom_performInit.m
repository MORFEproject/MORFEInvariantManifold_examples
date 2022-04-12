function [x0, v0] = HTHom_performInit(systemhandle, ~, param, activeParams, settings)
x0 = [];
v0 = [];

IP = settings.getSetting('IP');
s = IP.currentpointdata;
num = 1; %dummy
x = IP.source.x;
v = IP.source.v;

global HTHomds
HTHomds = settings.htds;
HTHomds.TestTolerance = settings.SParamTestTolerance;
gds = HTHomds.gds;
gds.T = HTHomds.T;
gds.eps1 = HTHomds.eps1;
ActiveSParams = find(gds.SParamsFree);
ActiveUParams = find(gds.UParamsFree);

sourcePoint = settings.coord;
gds.dim = length(sourcePoint);

    function reportError(message)
       assert(false, message); 
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START ORIGINAL CODE %%%%%%%%%%%%%%%%%%%%%%%%

number = 0; %number of sparams equal to zero
for i = 1:size(gds.SParams,1)
    if gds.SParams{i,2} == 0
        number = number+1;
    end
end

if number < (size(gds.SParams,1)-1)
    n_paru = size(ActiveUParams,2);
    n_pars = size(ActiveSParams,2);
    dim_npos = size(gds.UParams,1);
    
    for i = 1:n_pars
        if gds.SParams{ActiveSParams(i),2} == 0
            reportError('A zero SParam should not be denoted as free');
        end
    end
    
    if n_paru+n_pars ~= dim_npos+2
        reportError('The wrong number of free connection parameters is denoted');
    end
    
    n_par = size(activeParams,2);
    if n_par ~= 0
        reportError('The wrong number of free parameters is denoted');
    end
    
    if ~(gds.extravec(3)==1)
        reportError('eps1 must be free');
    end
    
    if gds.extravec(1)==1
        reportError('T must not be free');
    end
    
    if HTHomds.index == 0 %ConnA_ConnB
        HTHomds.index = 1;
        x = [x;sourcePoint];
        parU = vertcat(gds.UParams{:,2});
        %parU: alle UParams onder elkaar
        parS = vertcat(gds.SParams{:,2});
        %parS: ALLE SParams onder elkaar
        [x0,v0] = init_HTHom_HTHom(systemhandle,x,v,s, param, activeParams,parU,ActiveUParams,parS,ActiveSParams,settings.ntst,settings.ncol,gds.T,gds.eps1,settings.eps1tol);
    elseif HTHomds.index == 1 %ConnB_ConnB
        x0 = x(:,s(num).index);%enkel die éne kolom is belangrijk!
        x0 = [x0(1:(HTHomds.ntst*HTHomds.ncol+1)*length(sourcePoint));sourcePoint];
        v0 = v(:,s(num).index);
        parU = vertcat(gds.UParams{:,2});
        parS = vertcat(gds.SParams{:,2});
        [x0,v0] = init_HTHom_HTHom(systemhandle, x0, v0, s(num),HTHomds.P0,activeParams, parU, ActiveUParams, parS, ActiveSParams, settings.ntst,settings.ncol,gds.T,gds.eps1,settings.eps1tol);
    end
    
elseif number == (size(gds.SParams,1)-1)
    
    if ~(size(ActiveSParams,2) == 1)
        reportError('The wrong number of free connection parameters is denoted');
    end
    
    if gds.SParams{ActiveSParams(1),2} == 0
        reportError('A zero SParam should not be denoted as free');
    end
    
    if ~(size(ActiveUParams,2) == 0)
        reportError('The wrong number of free unstable parameters is denoted');
    end
    
    
    n_par = size(activeParams,2);
    if n_par ~= 1
        reportError('1 free system parameter is needed');
    end
    
    if ~(gds.extravec(3)==1)
        reportError('eps1 must be free');
    end
    
    if gds.extravec(1)==1
        reportError('T must not be free');
    end
    
    if HTHomds.index == 0 %ConnA_ConnC
        HTHomds.index = 2;
        parS = vertcat(gds.SParams{:,2});
        x = [x;sourcePoint];
        [x0,v0] = init_HTHom_HTHom(systemhandle, x, v, s, param, activeParams, [], [], parS, ActiveSParams, settings.ntst, settings.ncol,HTHomds.T,HTHomds.eps1,settings.eps1tol);
    elseif HTHomds.index == 1 %ConnB_ConnC
        HTHomds.index = 2;
        x0 = x(:,s(num).index);%enkel die éne kolom is belangrijk!
        v0 = v(:,s(num).index);
        x0 = [x0(1:(HTHomds.ntst*HTHomds.ncol+1)*length(sourcePoint));sourcePoint];%opgelet! sourcePoint ipv HTHomds.x0
        parS = vertcat(gds.SParams{:,2});
        [x0,v0] = init_HTHom_HTHom(systemhandle, x0, v0, s(num), param, activeParams, [], [], parS, ActiveSParams, settings.ntst, settings.ncol,gds.T,gds.eps1,settings.eps1tol);
    else %HTHomds.index == 2
        x0 = x(:,s(num).index);%enkel die éne kolom is belangrijk!
        v0 = v(:,s(num).index);
        x0 = [x0(1:(HTHomds.ntst*HTHomds.ncol+1)*length(sourcePoint));sourcePoint];
        parS = vertcat(gds.SParams{:,2});
        [x0,v0] = init_HTHom_HTHom(systemhandle, x0, v0, s(num), param, activeParams, [], [], parS, ActiveSParams, settings.ntst, settings.ncol,gds.T,gds.eps1,settings.eps1tol);
    end
    
else %ConnC_ConnD of ConnD_ConnD
    
    if ~(size(ActiveUParams,2) == 0)
        reportError('The wrong number of free unstable parameters is denoted');
    end
    
    if ~(size(ActiveSParams,2) == 0)
        reportError('The wrong number of free stable parameters is denoted');
    end
    
    n_par = size(activeParams,2);
    if n_par ~= 1
        reportError('1 free system parameter is needed');
    end
    
    if ~(gds.extravec(1)==1)
        reportError('T must be free');
    end
    if ~(gds.extravec(3)==1)
        reportError('eps1 must be free');
    end
    
    
    if HTHomds.index == 2 %ConnC_ConnD
        HTHomds.index = 3;
        x0 = x(:,s(num).index);
        v0 = v(:,s(num).index);
        x0 = [x0(1:(HTHomds.ntst*HTHomds.ncol+1)*length(sourcePoint));sourcePoint];
        gds.T = gds.period/2;   %NN: period is non-important?
        [x0,v0] = init_HTHom_HTHom(systemhandle, x0, v0, s(num), param, activeParams,[],[],[],[],settings.ntst,settings.ncol,gds.T,gds.eps1,settings.eps1tol);
    else %ConnD_ConnD
        x0 = x(:,s(num).index);
        v0 = v(:,s(num).index);
        x0 = [x0(1:(HTHomds.ntst*HTHomds.ncol+1)*length(sourcePoint));sourcePoint];
        gds.T = gds.period/2; %NN: period is non-important?
        [x0,v0] = init_HTHom_HTHom(systemhandle, x0, v0, s(num), param, activeParams,[],[],[],[],settings.ntst,settings.ncol,gds.T,gds.eps1,settings.eps1tol);
    end
    
end

end
