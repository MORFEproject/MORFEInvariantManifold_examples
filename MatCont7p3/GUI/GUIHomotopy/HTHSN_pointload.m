function settings = HTHSN_pointload(solution, settings, index)


htdssetting = settings.getSetting('htds');
gds = htdssetting.htds.gds;
HTHSNds = solution.globals.HTHSNds;

global cds;
cds = solution.globals.cds;   %restore cds to global variables for 'cjac';
x = solution.x;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

if HTHSNds.index == 1

    gds.T = HTHSNds.T;
    gds.eps0 = HTHSNds.eps0;
    gds.x0 = HTHSNds.x0;
    
    gds.UParams = []; 
    for i = 1:size(HTHSNds.UParams,1)
        gds.UParams{i,1} = strcat('UParam',num2str(i));        
        gds.UParams{i,2} = HTHSNds.UParams(i);
    end
    for i = 1:size(HTHSNds.ActiveUParams,2)
        gds.UParams{HTHSNds.ActiveUParams(i),2} = x(HTHSNds.ncoords+i,index);        
    end
    gds.SParams = []; 
    for i = 1:size(HTHSNds.SParams,1)
        gds.SParams{i,1} = strcat('SParam',num2str(i));        
        gds.SParams{i,2} = HTHSNds.SParams(i);
    end
    for i = 1:size(HTHSNds.ActiveSParams,2)
        gds.SParams{HTHSNds.ActiveSParams(i),2} = x(HTHSNds.ncoords+size(HTHSNds.ActiveUParams,2)+i,index);               
    end

    gds.eps1 = x(end,index);    
    HTHSNds.eps1 = gds.eps1; %NN
            
else
    %NN: checked.        
    ndim = size(HTHSNds.P0,1);
    for i = 1:ndim
        gds.parameters{i,2} = HTHSNds.P0(i,1);   
    end
    
    newind = HTHSNds.ncoords+1;
    
    gds.T = x(newind,index);
    newind = newind+1;
    gds.eps1 = x(newind,index);
    gds.period = gds.T*2;
    gds.x0 = HTHSNds.x0;
    gds.eps0 = HTHSNds.eps0;
    
    p= vertcat(gds.parameters{:,2});
    A = cjac(HTHSNds.func,HTHSNds.Jacobian,gds.x0,num2cell(p),[]);
    D = eig(A);

    [V,i] = min(abs(D));
    B = D([1:i-1 i+1:end]);
    nneg = sum(real(B) < 0);
    if (nneg == HTHSNds.nphase-1)
        if min(abs(real(B))) < 1e-10
            nneg = nneg -1;    
        end
    end
    if (nneg == 0)
        if min(abs(real(B))) < 1e-10
            nneg = nneg +1;
        end
    end

    npos = HTHSNds.nphase-nneg-1;   
    
    if nneg == 0
        QS = eye(HTHSNds.nphase);
    else
        [evc0,evl0] = eig(A);
        [K,i] = min(abs(diag(evl0)));
        evc = evc0(:,[1:i-1 i+1:end]);%eigenvectoren zonder eigenvector horende bij 0
        evl = evl0([1:i-1 i+1:end],[1:i-1 i+1:end]);%eigenwaarden zonder 0
    
        val = 0;
        pos = 0;
        k = 1;
        for i=1:size(A,1)-1
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
    
    [V,D] = eig(A);
    [Y,i] = min(abs(diag(D)));
    gds.SParams = [];
    if npos~=0
        Qtot = [QS(:,1:nneg) V(:,i)];
        [Qn,Rn] = qr(Qtot);
        sparams = zeros(1,npos);
        for i=1:npos        
            sparams(i) = 1/gds.eps1*(x(HTHSNds.ncoords-HTHSNds.nphase+1:HTHSNds.ncoords,index)-gds.x0)'*QS(:,nneg+1+i);
        end          

        for i = 1:length(sparams)
            gds.SParams{i,1} = strcat('SParam',num2str(i));
            gds.SParams{i,2} = sparams(i);
        end
    end

    if npos == 0    
        QU = eye(HTHSNds.nphase);                                                      
    else
        [evc0,evl0] = eig(A);
        [K,i] = min(abs(diag(evl0)));
        evc = evc0(:,[1:i-1 i+1:end]);%eigenvectoren zonder eigenvector horende bij 0
        evl = evl0([1:i-1 i+1:end],[1:i-1 i+1:end]);%eigenwaarden zonder 0
                       
        val = 0;                            
        pos = 0;                            
        k = 1;                            
        for i=1:size(A,1)-1                            
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

    [V,D] = eig(A);
    [Y,i] = min(abs(diag(D)));
    gds.UParams = [];
    Qtot = [QU(:,1:npos) V(:,i)];    
    [Qn,Rn] = qr(Qtot);    
    uparams = zeros(1,npos+1);    
    for i = 1:(npos+1)    
        uparams(i) = 1/gds.eps0*(x(1:HTHSNds.nphase,index)-gds.x0)'*Qn(:,i);        
    end
    
    som = sum(uparams.^2);            
    if abs(som-1)>1e-3    
        assert(false, 'It is not possible to start from this point.')                           
    end

    for i=1:length(uparams)    
        gds.UParams{i,1} = strcat('UParam',num2str(i));        
        gds.UParams{i,2} = uparams(i);        
    end

    if ~isfield(HTHSNds,'TestTolerance') || isempty(HTHSNds.TestTolerance)
        HTHSNds.TestTolerance = 1e-5;
    end

    if npos>1
        vec = [];
        vec = abs(vertcat(gds.SParams{:,2})) < HTHSNds.TestTolerance;    
        aantal = sum(vec);
        if aantal > 0
            HTHSNds.index = 1;        
        end
    end             
    HTHSNds.eps1 = gds.eps1;  %NN
    HTHSNds.T = gds.T; %NN   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NN
cds = [];
try gds = rmfield(gds, 'eps0'); catch; end
try gds = rmfield(gds, 'eps1'); catch; end
try gds = rmfield(gds, 'x0'); catch; end
try gds = rmfield(gds, 'parameters'); catch; end
try gds = rmfield(gds, 'coordinates'); catch; end
try gds = rmfield(gds, 'T'); catch; end

HTHSNds.gds = gds;
htdssetting.setValue(HTHSNds);
htdssetting.cleanupSParam(settings);
