function settings = HTHom_pointload(solution, settings, index)




htdssetting = settings.getSetting('htds');
gds = htdssetting.htds.gds;
HTHomds = solution.globals.HTHomds;

global cds;
cds = solution.globals.cds;  %restore cds to global variables for 'cjac';
x = solution.x;

if HTHomds.index == 1
    
    gds.T = HTHomds.T;
    gds.eps0 = HTHomds.eps0;
    gds.x0 = HTHomds.x0;
    
    for i = 1:size(HTHomds.ActiveUParams,2)
        gds.UParams{HTHomds.ActiveUParams(i),2} = x(HTHomds.nphase*(HTHomds.ntst*HTHomds.ncol+1)+i,index);        
    end
    for i = 1:size(HTHomds.ActiveSParams,2)
        gds.SParams{HTHomds.ActiveSParams(i),2} = x(HTHomds.nphase*(HTHomds.ntst*HTHomds.ncol+1)+size(HTHomds.ActiveUParams,2)+i,index);               
    end

    gds.eps1 = x(end,index);  
    HTHomds.eps1 = gds.eps1; %NN
    
elseif HTHomds.index == 2 
    
 
    gds.x0 = x(HTHomds.nphase*HTHomds.tps+1:HTHomds.nphase*(HTHomds.tps+1),index);
    ndim = size(HTHomds.P0,1);
    for i = 1:ndim
        gds.parameters{i,2} = HTHomds.P0(i,1);    
    end
    j = HTHomds.ActiveParams(1);
    gds.parameters{j,2} = x(HTHomds.nphase*(HTHomds.tps+1)+1,index);   
    
    gds.T = HTHomds.T;              

    p = vertcat(gds.parameters{:,2});
    A = cjac(HTHomds.func,HTHomds.Jacobian,gds.x0,num2cell(p),HTHomds.ActiveParams);
    D = eig(A);
    nneg = sum(real(D) < 0);
    
    % If one eigenvalue is (practically) zero, and the one of the subspaces has
    % zero dimension, change this dimension with 1.
    if (nneg == HTHomds.nphase)
        if min(abs(real(D))) < 1e-10
            nneg = nneg -1;
        end
    end
    if (nneg == 0)
        if min(abs(real(D))) < 1e-10
            nneg = nneg +1;
        end
    end
    npos = HTHomds.nphase-nneg;
    
    if nneg == 0
        QS = eye(HTHomds.nphase);
    else
        [evc,evl] = eig(A);
    
        val = 0;
        pos = 0;
        k = 1;
        for i=1:size(A,1)
            if real(evl(i,i)) < 0
                val(k) = real(evl(i,i));
                pos(k) = i; %geeft je de posities waar de gezochte eigenwaarden en eigenvectoren staan
                k = k+1;
            end
        end
        evcs = evc(:,pos); %de eigenvectoren die meetellen voor de negatieve
        evls = 0;
        for t = 1:size(val,2)
            evls(t) = evl(pos(t),pos(t));
        end

        [a,b] = sort(val); % dan is a(i) = val(b(i))
        VU = zeros(size(A,1),size(val,2));
        sortedevls = zeros(size(pos,2),1);
        for l = 1:size(val,2)
            VU(:,l) = evcs(:,b(l));
            sortedevls(l) = evls(b(l)); %de eigenwaarden gesorteerd van klein naar groot (enkel negatief reeel deel)
        end   

        VU_1 = zeros(size(A,1),size(val,2));
        sortedevls_1 = zeros(size(pos,2),1);
        for f = 1:size(val,2)
            VU_1(:,f) = VU(:,end-f+1);
            sortedevls_1(f) = sortedevls(end-f+1);
        end
        sortedevls = sortedevls_1;

        B = VU_1;
        k = 1;
        while k <= size(B,2)
            ok = 1;
            init = 1;
            while ok == 1 && init <= size(B,1)
                if imag(B(init,k))                     
                    tmp = B(:,k);
                    B(:,k) = real(tmp);
                    B(:,k+1) = imag(tmp);
                    ok = 0;
                end
                init = init+1;
            end
            if ok == 1
               k = k+1; 
            else
                k = k+2;
            end            
        end

        % Compute orthonormal basis for the eigenspace
        [QS,RS] = qr(B);
    end
        
    gds.SParams = [];    
    gds.eps1 = x(HTHomds.nphase*(HTHomds.tps+1)+3,index);
    if npos ~= 0
        sparams = zeros(1,npos);
        for i=1:npos        
            sparams(i) = 1/gds.eps1*(x(HTHomds.ncoords-HTHomds.nphase+1:HTHomds.ncoords,index)-gds.x0)'*QS(:,nneg+i);
        end          
        for i = 1:length(sparams)
            gds.SParams{i,1} = strcat('SParam',num2str(i));
            gds.SParams{i,2} = sparams(i);
        end
    end

    gds.eps1 = x(HTHomds.nphase*(HTHomds.tps+1)+3,index);
    gds.eps0 = HTHomds.eps0;
        
    if npos == 0    
        QU = eye(HTHomds.nphase);            
    else
        [evc,evl] = eig(A);
                       
        val = 0;        
        pos = 0;        
        k = 1;        
        for i=1:size(A,1)        
            if real(evl(i,i)) > 0            
                val(k) = real(evl(i,i));                
                pos(k) = i; %geeft je de posities waar de gezochte eigenwaarden en eigenvectoren staan                
                k = k+1;                
            end
        end
        evcs = evc(:,pos); %de eigenvectoren die meetellen        
        evls = 0;        
        for t = 1:size(val,2)        
            evls(t) = evl(pos(t),pos(t));            
        end

        [a,b] = sort(val); % dan is a(i) = val(b(i))        
        VU = zeros(size(A,1),size(val,2));        
        sortedevls = zeros(size(pos,2),1);        
        for l = 1:size(val,2)        
            VU(:,l) = evcs(:,b(l));            
            sortedevls(l) = evls(b(l)); %de eigenwaarden gesorteerd van klein naar groot (enkel positief reeel deel)            
        end

        B = VU;        
        k = 1;        
        while k <= size(B,2)        
            ok = 1;            
            init = 1;           
            while ok == 1 && init <= size(B,1)            
                if imag(B(init,k))                                     
                    tmp = B(:,k);                    
                    B(:,k) = real(tmp);                    
                    B(:,k+1) = imag(tmp);                    
                    ok = 0;                    
                end
                init = init+1;                
            end
            if ok == 1            
                k = k+1;                 
            else
                k = k+2;                
            end
        end

        % Compute orthonormal basis for the eigenspace        
        [QU,RU] = qr(B);        
    end
    
    gds.UParams = [];
    if npos ~= 0
        uparams = zeros(1,npos);
        for i = 1:npos
            uparams(i) = 1/gds.eps0*(x(1:HTHomds.nphase,index)-gds.x0)'*QU(:,i);
        end
    
        som = sum(uparams.^2);        
        if abs(som-1)>1e-3
            assert(false, 'It is not possible to start from this point.')                   
        end
                
        for i=1:length(uparams)
            gds.UParams{i,1} = strcat('UParam',num2str(i));
            gds.UParams{i,2} = uparams(i);
        end
    end                       


    if npos>1
        vec = [];
        vec = abs(vertcat(gds.SParams{:,2})) < HTHomds.TestTolerance;    
        aantal = sum(vec);
        if aantal < (npos-1) 
            HTHomds.index = 1;        
        end
    end
            
    HTHomds.eps1 = gds.eps1;  %NN
    HTHomds.T = gds.T; %NN
else


    gds.x0 = x(HTHomds.nphase*HTHomds.tps+1:HTHomds.nphase*(HTHomds.tps+1),index);
    
    ndim = size(HTHomds.P0,1);
    for i = 1:ndim
        gds.parameters{i,2} = HTHomds.P0(i,1);   
    end
    len = length(HTHomds.ActiveParams);
    for i = 1:len
        j = HTHomds.ActiveParams(1,i);    
        gds.parameters{j,2} = x(HTHomds.nphase*(HTHomds.tps+1)+i,index);    
    end
    newind = HTHomds.nphase*(HTHomds.tps+1)+len+1;

    gds.T = x(newind,index);
    newind = newind+1;
    gds.eps1 = x(newind,index);
    gds.period = gds.T*2;
    gds.eps0 = HTHomds.eps0;   
    
    p= vertcat(gds.parameters{:,2});
    A = cjac(HTHomds.func,HTHomds.Jacobian,gds.x0,num2cell(p),HTHomds.ActiveParams);
    D = eig(A);
    nneg = sum(real(D) < 0);
    
    % If one eigenvalue is (practically) zero, and the one of the subspaces has
    % zero dimension, change this dimension with 1.
    if (nneg == HTHomds.nphase)
        if min(abs(real(D))) < 1e-10
            nneg = nneg -1;
        end
    end
    if (nneg == 0)
        if min(abs(real(D))) < 1e-10
            nneg = nneg +1;
        end
    end
    npos = HTHomds.nphase-nneg;
    
    if nneg == 0
        QS = [];
    else
        [evc,evl] = eig(A);
    
        val = 0;
        pos = 0;
        k = 1;
        for i=1:size(A,1)
            if real(evl(i,i)) < 0
                val(k) = real(evl(i,i));
                pos(k) = i; %geeft je de posities waar de gezochte eigenwaarden en eigenvectoren staan
                k = k+1;
            end
        end
        evcs = evc(:,pos); %de eigenvectoren die meetellen voor de negatieve
        evls = 0;
        for t = 1:size(val,2)
            evls(t) = evl(pos(t),pos(t));
        end

        [a,b] = sort(val); % dan is a(i) = val(b(i))
        VU = zeros(size(A,1),size(val,2));
        sortedevls = zeros(size(pos,2),1);
        for l = 1:size(val,2)
            VU(:,l) = evcs(:,b(l));
            sortedevls(l) = evls(b(l)); %de eigenwaarden gesorteerd van klein naar groot (enkel negatief reeel deel)
        end   

        VU_1 = zeros(size(A,1),size(val,2));
        sortedevls_1 = zeros(size(pos,2),1);
        for f = 1:size(val,2)
            VU_1(:,f) = VU(:,end-f+1);
            sortedevls_1(f) = sortedevls(end-f+1);
        end
        sortedevls = sortedevls_1;

        B = VU_1;
        k = 1;
        while k <= size(B,2)
            ok = 1;
            init = 1;
            while ok == 1 && init <= size(B,1)
                if imag(B(init,k))                     
                    tmp = B(:,k);
                    B(:,k) = real(tmp);
                    B(:,k+1) = imag(tmp);
                    ok = 0;
                end
                init = init+1;
            end
            if ok == 1
               k = k+1; 
            else
                k = k+2;
            end            
        end

        % Compute orthonormal basis for the eigenspace
        [QS,RS] = qr(B);
    end
    
    
    gds.SParams = [];
    if npos~=0
        sparams = zeros(1,npos);
        for i=1:npos        
            sparams(i) = 1/gds.eps1*(x(HTHomds.ncoords-HTHomds.nphase+1:HTHomds.ncoords,index)-gds.x0)'*QS(:,nneg+i);
        end          

        for i = 1:length(sparams)
            gds.SParams{i,1} = strcat('SParam',num2str(i));
            gds.SParams{i,2} = sparams(i);
        end
    end

    if npos == 0    
        QU = [];                            
        sortedevls = [];                            
    else
        [evc,evl] = eig(A);
                       
        val = 0;                            
        pos = 0;                            
        k = 1;                            
        for i=1:size(A,1)                            
            if real(evl(i,i)) > 0                                        
                val(k) = real(evl(i,i));                                                    
                pos(k) = i; %geeft je de posities waar de gezochte eigenwaarden en eigenvectoren staan                                                    
                k = k+1;                                                    
            end
        end
        evcs = evc(:,pos); %de eigenvectoren die meetellen        
        evls = 0;                           
        for t = 1:size(val,2)                            
            evls(t) = evl(pos(t),pos(t));                                        
        end

        [a,b] = sort(val); % dan is a(i) = val(b(i))        
        VU = zeros(size(A,1),size(val,2));                            
        sortedevls = zeros(size(pos,2),1);                           
        for l = 1:size(val,2)                            
            VU(:,l) = evcs(:,b(l));                                        
            sortedevls(l) = evls(b(l)); %de eigenwaarden gesorteerd van klein naar groot (enkel positief reeel deel)                                        
        end

        B = VU;        
        k = 1;                            
        while k <= size(B,2)                            
            ok = 1;                                        
            init = 1;                                       
            while ok == 1 && init <= size(B,1)                                       
                if imag(B(init,k))                                                                         
                    tmp = B(:,k);                                                                
                    B(:,k) = real(tmp);                                                                
                    B(:,k+1) = imag(tmp);                                                                
                    ok = 0;                                                                
                end
                init = init+1;               
            end
            if ok == 1            
                k = k+1;                                                     
            else
                k = k+2;                
            end
        end
        % Compute orthonormal basis for the eigenspace        
        [QU,RU] = qr(B);                            
    end

    
    gds.UParams = [];
    if npos ~= 0
        uparams = zeros(1,npos);
        for i = 1:npos
            uparams(i) = 1/gds.eps0*(x(1:HTHomds.nphase,index)-gds.x0)'*QU(:,i);
        end
        som = sum(uparams.^2);      
        
        if abs(som-1)>1e-3
            assert(false, 'It is not possible to start from this point.')                   
        end
                
        for i=1:length(uparams)
            gds.UParams{i,1} = strcat('UParam',num2str(i));
            gds.UParams{i,2} = uparams(i);
        end
    end                       

    
    if npos>1
        vec = [];
        vec = abs(vertcat(gds.SParams{:,2})) < HTHomds.TestTolerance;    
        aantal = sum(vec);
        if aantal < (npos-1) 
            HTHomds.index = 1;        
        elseif ~(aantal==(HTHomds.npos))
            HTHomds.index = 2;
        end
    end             
    HTHomds.eps1 = gds.eps1;  %NN
    HTHomds.T = gds.T; %NN
end

%NN
cds = [];
try gds = rmfield(gds, 'eps0'); catch; end
try gds = rmfield(gds, 'eps1'); catch; end
try gds = rmfield(gds, 'x0'); catch; end
try gds = rmfield(gds, 'parameters'); catch; end
try gds = rmfield(gds, 'coordinates'); catch; end
try gds = rmfield(gds, 'T'); catch; end

HTHomds.gds = gds;
htdssetting.setValue(HTHomds);
htdssetting.cleanupSParam(settings);
