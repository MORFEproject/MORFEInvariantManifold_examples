function settings = HTHet_pointload(solution, settings, index)


htdssetting = settings.getSetting('htds');
gds = htdssetting.htds.gds;
HTHetds = solution.globals.HTHetds;

global cds;
cds = solution.globals.cds;   %restore cds to global variables for 'cjac';
x = solution.x;

dim = length(HTHetds.x0);
ndim = length(HTHetds.P0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if HTHetds.index == 1
    for i = 1:dim
        gds.coordinates{i,2} = x(i,index);
    end
    for i = 1:ndim
        gds.parameters{i,2} = HTHetds.P0(i,1);
    end
    
    gds.T = HTHetds.T;
    gds.eps0 = HTHetds.eps0;
    gds.x0 = HTHetds.x0;
    gds.x1 = HTHetds.x1;
    
    gds.UParams = [];
    for i = 1:length(HTHetds.UParams)
        gds.UParams{i,1} = strcat('UParam',num2str(i));
        gds.UParams{i,2} = HTHetds.UParams(i);
    end
    for i = 1:size(HTHetds.ActiveUParams,2)
        gds.UParams{HTHetds.ActiveUParams(i),2} = x(HTHetds.ncoords+i,index);
    end
    
    gds.SParams = [];
    for i = 1:length(HTHetds.SParams)
        gds.SParams{i,1} = strcat('SParam',num2str(i));
        gds.SParams{i,2} = HTHetds.SParams(i);
    end
    for i = 1:size(HTHetds.ActiveSParams,2)
        gds.SParams{HTHetds.ActiveSParams(i),2} = x(HTHetds.ncoords+size(HTHetds.ActiveUParams,2)+i,index);
    end
    
    gds.eps1 = x(end,index);
    HTHetds.T = gds.T;
    HTHetds.eps1 = gds.eps1;
    HTHetds.eps0 = gds.eps0;   
elseif HTHetds.index == 2
    
    for i = 1:dim
        gds.coordinates{i,2} = x(i,index);
    end
    
    gds.x0 = x(HTHetds.ncoords+1:HTHetds.ncoords+HTHetds.nphase,index);
    gds.x1 = x(HTHetds.ncoords+HTHetds.nphase+1:HTHetds.ncoords+2*HTHetds.nphase,index);
    ndim = size(HTHetds.P0,1);
    for i = 1:ndim
        gds.parameters{i,2} = HTHetds.P0(i,1);
    end
    for i = 1:length(HTHetds.ActiveParams)
        gds.parameters{HTHetds.ActiveParams(i),2} = x(HTHetds.PeriodIdx+i,index);
    end
    
    gds.T = HTHetds.T;
    gds.eps1 = x(HTHetds.PeriodIdx+length(HTHetds.ActiveParams)+length(HTHetds.ActiveSParams)+1,index);
    gds.eps0 = HTHetds.eps0;
    gds.period = gds.T*2;
    
    p = vertcat(gds.parameters{:,2});
    A = cjac(HTHetds.func,HTHetds.Jacobian,gds.x1,num2cell(p),HTHetds.ActiveParams);
    D = eig(A);
    nneg = sum(real(D) < 0);
    
    % If one eigenvalue is (practically) zero, and the one of the subspaces has
    % zero dimension, change this dimension with 1.
    if (nneg == HTHetds.nphase)
        if min(abs(real(D))) < 1e-10
            nneg = nneg -1;
        end
    end
    if (nneg == 0)
        if min(abs(real(D))) < 1e-10
            nneg = nneg +1;
        end
    end
    
    if nneg == 0
        QS = eye(HTHetds.nphase);
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
    
    %    QbS1 = QS(:,1:nneg);
    %    QbS2 = QS(:,nneg+1:end);
    %    YS = x(end-nneg*(HTHetds.nphase-nneg)+1:end,index);
    
    %    [Q1S,S1,R1] = svd(QbS1 + QbS2*YS);
    %    Q1S = Q1S*diag(sign(diag(HTHetds.oldStableQ'*Q1S)));
    
    %    reshape(YS,HTHetds.nphase-nneg,nneg);
    %    Q1S = QS * [-YS'; eye(size(YS,1))];
    
    gds.SParams = [];
    if (HTHetds.nphase-nneg) ~= 0
        sparams = zeros(1,HTHetds.nphase-nneg);
        for i=1:(HTHetds.nphase-nneg)
            sparams(i) = 1/gds.eps1*(x(HTHetds.ncoords-HTHetds.nphase+1:HTHetds.ncoords,index)-gds.x1)'*QS(:,nneg+i);
        end
        
        for i = 1:length(sparams)
            gds.SParams{i,1} = strcat('SParam',num2str(i));
            gds.SParams{i,2} = sparams(i);
        end
    end
    
    p = vertcat(gds.parameters{:,2});
    A = cjac(HTHetds.func,HTHetds.Jacobian,gds.x0,num2cell(p),HTHetds.ActiveParams);
    D = eig(A);
    npos = sum(real(D) > 0);
    
    % If one eigenvalue is (practically) zero, and the one of the subspaces has
    % zero dimension, change this dimension with 1.
    if (npos == HTHetds.nphase)
        if min(abs(real(D))) < 1e-10
            npos = npos -1;
        end
    end
    if (npos == 0)
        if min(abs(real(D))) < 1e-10
            npos = npos +1;
        end
    end
    
    if npos == 0
        QU = [];
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
            uparams(i) = 1/gds.eps0*(x(1:HTHetds.nphase,index)-gds.x0)'*QU(:,i);
        end
        
        som = sum(uparams.^2);
        if abs(som-1)>1e-3
            error('It is not possible to start from this point.')
        end
        
        for i=1:length(uparams)
            gds.UParams{i,1} = strcat('UParam',num2str(i));
            gds.UParams{i,2} = uparams(i);
        end
    end
    
    if ~isfield(HTHetds,'TestTolerance') || isempty(HTHetds.TestTolerance)
        HTHetds.TestTolerance = 1e-5;
    end
    
    if npos>1
        vec = [];
        vec = abs(vertcat(gds.SParams{:,2})) < HTHetds.TestTolerance;
        aantal = sum(vec);
        if aantal < (npos-1)
            HTHetds.index = 1;
        end
    end
    
    HTHetds.T = gds.T;
    HTHetds.eps1 = gds.eps1;
    HTHetds.eps0 = gds.eps0;
else
    
    for i = 1:dim
        gds.coordinates{i,2} = x(i,index);
    end
    gds.x0 = x(HTHetds.ncoords+1:HTHetds.ncoords+HTHetds.nphase,index);
    gds.x1 = x(HTHetds.ncoords+HTHetds.nphase+1:HTHetds.ncoords+2*HTHetds.nphase,index);
    
    ndim = size(HTHetds.P0,1);
    for i = 1:ndim
        gds.parameters{i,2} = HTHetds.P0(i,1);
    end
    for i = 1:length(HTHetds.ActiveParams)
        gds.parameters{HTHetds.ActiveParams(i),2} = x(HTHetds.PeriodIdx+i,index);
    end
    newind = HTHetds.PeriodIdx+length(HTHetds.ActiveParams)+1;
    
    gds.T = x(newind,index);
    newind = newind+1;
    gds.eps1 = x(newind,index);
    gds.period = gds.T*2;
    
    p= vertcat(gds.parameters{:,2});
    A = cjac(HTHetds.func,HTHetds.Jacobian,gds.x1,num2cell(p),HTHetds.ActiveParams);
    D = eig(A);
    nneg = sum(real(D) < 0);
    
    % If one eigenvalue is (practically) zero, and the one of the subspaces has
    % zero dimension, change this dimension with 1.
    if (nneg == HTHetds.nphase)
        if min(abs(real(D))) < 1e-10
            nneg = nneg -1;
        end
    end
    if (nneg == 0)
        if min(abs(real(D))) < 1e-10
            nneg = nneg +1;
        end
    end
    
    
    if (HTHetds.nphase-nneg) == 0
        QS = eye(HTHetds.nphase);
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
    if (HTHetds.nphase-nneg) ~= 0
        sparams = zeros(1,HTHetds.nphase-nneg);
        for i=1:(HTHetds.nphase-nneg)
            sparams(i) = 1/gds.eps1*(x(HTHetds.ncoords-HTHetds.nphase+1:HTHetds.ncoords,index)-gds.x1)'*QS(:,nneg+i);
        end
        
        for i = 1:length(sparams)
            gds.SParams{i,1} = strcat('SParam',num2str(i));
            gds.SParams{i,2} = sparams(i);
        end
    end
    
    
    p = vertcat(gds.parameters{:,2});
    A = cjac(HTHetds.func,HTHetds.Jacobian,gds.x0,num2cell(p),HTHetds.ActiveParams);
    D = eig(A);
    npos = sum(real(D) > 0);
    
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
            uparams(i) = 1/HTHetds.eps0*(x(1:HTHetds.nphase,index)-gds.x0)'*QU(:,i);
        end
        
        som = sum(uparams.^2);
        if abs(som-1)>1e-3
            error('It is not possible to start from this point.')
        end
        
        for i=1:length(uparams)
            gds.UParams{i,1} = strcat('UParam',num2str(i));
            gds.UParams{i,2} = uparams(i);
        end
    end
    
    if ~isfield(HTHetds,'TestTolerance') || isempty(HTHetds.TestTolerance)
        HTHetds.TestTolerance = 1e-5;
    end
    
    vec = [];
    vec = abs(vertcat(gds.SParams{:,2})) < HTHetds.TestTolerance;
    aantal = sum(vec);
    if npos>1
        if aantal < (npos-1)
            HTHetds.index = 1;
        end
    end
    if (npos>1) && (aantal>=npos-1) && ~(aantal==(HTHetds.nphase-nneg))
        HTHetds.index = 2;
    end
    HTHetds.T = gds.T;
    HTHetds.eps1 = gds.eps1;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NN
settings.coord.set(gds.x0);
settings.coord_target.set(gds.x1);

cds = [];
try gds = rmfield(gds, 'eps0'); catch; end
try gds = rmfield(gds, 'eps1'); catch; end
try gds = rmfield(gds, 'x0'); catch; end
try gds = rmfield(gds, 'x1'); catch; end
try gds = rmfield(gds, 'parameters'); catch; end
try gds = rmfield(gds, 'coordinates'); catch; end
try gds = rmfield(gds, 'T'); catch; end

HTHetds.gds = gds;
htdssetting.setValue(HTHetds);
htdssetting.cleanupSParam(settings);
