function [x0, v0] = HTHSN_performInit(systemhandle, ~, param, activeParams, settings)
x0 = [];
v0 = [];

IP = settings.getSetting('IP');
s = IP.currentpointdata;
num = 1; %dummy
x = IP.source.x;
v = IP.source.v;

global HTHSNds
HTHSNds = settings.htds;
HTHSNds.TestTolerance = settings.SParamTestTolerance;
gds = HTHSNds.gds;
gds.T = HTHSNds.T;
gds.eps1 = HTHSNds.eps1;
gds.dim = length(settings.coord);
ActiveSParams = find(gds.SParamsFree);
ActiveUParams = find(gds.UParamsFree);

sourcePoint = settings.coord;


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

if (size(gds.SParams,1)-number) > 0
    n_paru = size(ActiveUParams,2);
    n_pars = size(ActiveSParams,2);
    dim_npos = size(gds.UParams,1)-1;
    
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
    
    if HTHSNds.index == 0
        HTHSNds.index = 1;
        parU = vertcat(gds.UParams{:,2});
        %parU: alle UParams onder elkaar
        parS = vertcat(gds.SParams{:,2});
        %parS: ALLE SParams onder elkaar
        x = [x;sourcePoint];
        par = param;
        [x0,v0] = init_HTHSN_HTHSN(systemhandle,x,v,s,par,parU,ActiveUParams,parS,ActiveSParams,settings.ntst,settings.ncol,gds.T,gds.eps1,settings.eps1tol);
    elseif HTHSNds.index == 1
        x0 = x(:,s(num).index);%enkel die éne kolom is belangrijk!
        x0 = [x0(1:(HTHSNds.ntst*HTHSNds.ncol+1)*gds.dim);sourcePoint];
        v0 = v(:,s(num).index);
        parU = vertcat(gds.UParams{:,2});
        parS = vertcat(gds.SParams{:,2});
        par = param;
        [x0,v0] = init_HTHSN_HTHSN(systemhandle, x0, v0, s(num),par,parU, ActiveUParams, parS, ActiveSParams, settings.ntst,settings.ncol,gds.T,gds.eps1,settings.eps1tol);
    end
    
else
    
    if ~(size(ActiveUParams,2) == 0)
        reportError('The wrong number of free unstable parameters is denoted');
    end
    
    if ~(size(ActiveSParams,2) == 0)
        reportError('The wrong number of free stable parameters is denoted');
    end
    
    n_par = size(activeParams,2);
    if n_par ~= 0
        reportError('0 free system parameters are needed');
    end
    
    if ~(gds.extravec(1)==1)
        reportError('T must be free');
    end
    if ~(gds.extravec(3)==1)
        reportError('eps1 must be free');
    end
    
    if HTHSNds.index == 0
        HTHSNds.index = 2;
        x = [x;sourcePoint];
        par = param;
        [x0,v0] = init_HTHSN_HTHSN(systemhandle,x,v,s,par,[],[],[],[],settings.ntst,settings.ncol,gds.T,gds.eps1,settings.eps1tol);
    elseif HTHSNds.index == 1
        HTHSNds.index = 2;
        x0 = x(:,s(num).index);
        v0 = v(:,s(num).index);
        x0 = [x0(1:(HTHSNds.ntst*HTHSNds.ncol+1)*gds.dim);sourcePoint];
        par = param;
        gds.T = gds.period/2;
        [x0,v0] = init_HTHSN_HTHSN(systemhandle, x0, v0, s(num), par,[],[],[],[],settings.ntst,settings.ncol,gds.T,gds.eps1,settings.eps1tol);
    else
        x0 = x(:,s(num).index);
        v0 = v(:,s(num).index);
        x0 = [x0(1:(HTHSNds.ntst*HTHSNds.ncol+1)*gds.dim);sourcePoint];
        par = param;
        gds.T = gds.period/2;
        [x0,v0] = init_HTHSN_HTHSN(systemhandle, x0, v0, s(num), par,[],[],[],[],settings.ntst,settings.ncol,gds.T,gds.eps1,settings.eps1tol);
    end
    
end



end