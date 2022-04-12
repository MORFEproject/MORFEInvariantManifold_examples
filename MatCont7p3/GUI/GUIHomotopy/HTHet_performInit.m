function [x0, v0] = HTHet_performInit(systemhandle, ~, param, activeParams, settings)
x0 = [];
v0 = [];

IP = settings.getSetting('IP');
s = IP.currentpointdata;
num = 1; %dummy
x = IP.source.x;
v = IP.source.v;

global HTHetds
HTHetds = settings.htds;
HTHetds.TestTolerance = settings.SParamTestTolerance;
gds = HTHetds.gds;
gds.T = HTHetds.T;
gds.eps1 = HTHetds.eps1;
gds.dim = length(settings.coord);
ActiveSParams = find(gds.SParamsFree);
ActiveUParams = find(gds.UParamsFree);

sourcePoint = settings.coord;
targetPoint = settings.coord_target;
gds.x0 = sourcePoint(:);
gds.x1 = targetPoint(:);


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

dim_npos = size(gds.UParams,1);
tmp = find(vertcat(gds.SParams{:,2})~=0);

if ~isempty(tmp)
    if ~(length(tmp)==length(ActiveSParams)) || ~isequal(tmp',ActiveSParams)  
        reportError('The wrong SParams are denoted as free');        
    end
elseif ~(length(tmp)==length(ActiveSParams))                
    reportError('The wrong SParams are denoted as free');        
end      

if number < (dim_npos-1)
    n_paru = size(ActiveUParams,2);
    n_pars = size(ActiveSParams,2);
    
    for i = 1:n_pars
        if gds.SParams{ActiveSParams(i),2} == 0                
            reportError('A zero SParam should not be denoted as free');
        end        
    end          
    
    n_par = size(activeParams,2);
    if n_par ~= 0
        reportError('0 parameters should be denoted as free');
    end
                               
    if ~(gds.extravec(3)==1)
        reportError('eps1 must be free');
    end
    
    if gds.extravec(1)==1
        reportError('T must not be free');
    end        
  
    if HTHetds.index == 0 %ConnA_ConnB      
        HTHetds.index = 1;
        parU = vertcat(gds.UParams{:,2});            
        %parU: alle UParams onder elkaar            
        parS = vertcat(gds.SParams{:,2});                            
        %parS: ALLE SParams onder elkaar     
        x = [x;gds.x0;gds.x1];
        [x0,v0] = init_HTHet_HTHet(systemhandle,x,v,s,param, activeParams,parU,ActiveUParams,parS,ActiveSParams,settings.ntst,settings.ncol,gds.T,gds.eps1,settings.eps1tol);
    elseif HTHetds.index == 1 %ConnB_ConnB           
        x0 = x(:,s(num).index);%enkel die éne kolom is belangrijk!                                  
        x0 = [x0(1:(HTHetds.ntst*HTHetds.ncol+1)*gds.dim);gds.x0;gds.x1];
        v0 = v(:,s(num).index);            
        parU = vertcat(gds.UParams{:,2});            
        parS = vertcat(gds.SParams{:,2});              
        [x0,v0] = init_HTHet_HTHet(systemhandle, x0, v0, s(num),param,activeParams, parU, ActiveUParams, parS, ActiveSParams, settings.ntst,settings.ncol,gds.T,gds.eps1,settings.eps1tol);            
    end
    
elseif (size(gds.SParams,1)-number) > 0        
           
    for i = 1:size(ActiveSParams,2)
        if gds.SParams{ActiveSParams(i),2} == 0                
            reportError('A zero SParam should not be denoted as free');
        end
    end
    
    if ~(size(ActiveUParams,2) == 0)
        reportError('The wrong number of free unstable parameters is denoted');
    end                             
    
    if ~(gds.extravec(3)==1)
        reportError('eps1 must be free');
    end
    
    if gds.extravec(1)==1
        reportError('T must not be free');
    end
    
    if HTHetds.index == 0 %ConnA_ConnC            
        HTHetds.index = 2;
        parS = vertcat(gds.SParams{:,2});
        x = [x;gds.x0;gds.x1];
        [x0,v0] = init_HTHet_HTHet(systemhandle, x, v, s, param, activeParams, [], [], parS, ActiveSParams, settings.ntst, settings.ncol,gds.T,gds.eps1,settings.eps1tol);
    elseif HTHetds.index == 1 %ConnB_ConnC     
        HTHetds.index = 2;
        x0 = x(:,s(num).index);%enkel die éne kolom is belangrijk!     
        v0 = v(:,s(num).index);    
        x0 = [x0(1:(HTHetds.ntst*HTHetds.ncol+1)*gds.dim);gds.x0;gds.x1];%opgelet! gds.x0 ipv HTHetds.x0
        par = param;          
        parS = vertcat(gds.SParams{:,2});           
        [x0,v0] = init_HTHet_HTHet(systemhandle, x0, v0, s(num), par, activeParams, [], [], parS, ActiveSParams, settings.ntst, settings.ncol,gds.T,gds.eps1,settings.eps1tol);      
    else %HTHetds.index == 2       
        x0 = x(:,s(num).index);%enkel die éne kolom is belangrijk!                    
        v0 = v(:,s(num).index);    
        x0 = [x0(1:(HTHetds.ntst*HTHetds.ncol+1)*gds.dim);gds.x0;gds.x1];
        par = param;          
        parS = vertcat(gds.SParams{:,2});
        [x0,v0] = init_HTHet_HTHet(systemhandle, x0, v0, s(num), par, activeParams, [], [], parS, ActiveSParams, settings.ntst, settings.ncol,gds.T,gds.eps1,settings.eps1tol);                                    
    end           
    
else %ConnC_ConnD of ConnD_ConnD
    
    if ~(size(ActiveUParams,2) == 0)
        reportError('0 unstable parameters should be denoted as free');
    end
    
    if ~(size(ActiveSParams,2) == 0)
        reportError('0 stable parameters should be denoted as free');
    end        
    
    if ~(gds.extravec(1)==1)
        reportError('T must be free');
    end       
    if ~(gds.extravec(3)==1)
        reportError('eps1 must be free');
    end
    

    if HTHetds.index == 1
        HTHetds.index = 3;
        x0 = x(:,s(num).index);            
        v0 = v(:,s(num).index);       
        x0 = [x0(1:(HTHetds.ntst*HTHetds.ncol+1)*gds.dim);gds.x0;gds.x1];
        par = param;            
        gds.T = gds.period/2;            
        [x0,v0] = init_HTHet_HTHet(systemhandle, x0, v0, s(num), par, activeParams,[],[],[],[],settings.ntst,settings.ncol,gds.T,gds.eps1,settings.eps1tol);            
    elseif HTHetds.index == 2 %ConnC_ConnD
        HTHetds.index = 3;
        x0 = x(:,s(num).index);            
        v0 = v(:,s(num).index);                
        x0 = [x0(1:(HTHetds.ntst*HTHetds.ncol+1)*gds.dim);gds.x0;gds.x1];
        par = param;            
        gds.T = gds.period/2;            
        [x0,v0] = init_HTHet_HTHet(systemhandle, x0, v0, s(num), par, activeParams,[],[],[],[],settings.ntst,settings.ncol,gds.T,gds.eps1,settings.eps1tol);            
    else %ConnD_ConnD
        x0 = x(:,s(num).index);            
        v0 = v(:,s(num).index);            
        x0 = [x0(1:(HTHetds.ntst*HTHetds.ncol+1)*gds.dim);gds.x0;gds.x1];
        par = param;            
        gds.T = gds.period/2;            
        [x0,v0] = init_HTHet_HTHet(systemhandle, x0, v0, s(num), par, activeParams,[],[],[],[],settings.ntst,settings.ncol,gds.T,gds.eps1,settings.eps1tol);
    end
    
end                       


end
