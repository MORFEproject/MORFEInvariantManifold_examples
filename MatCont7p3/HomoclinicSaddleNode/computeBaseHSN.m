function  [Q0,sortedevls] =  computeBaseHSN(A0,flag,dim)
%flag == 0: unstable, flag ==1: stable
global homds

if flag == 0
    if dim == 0
        Q0 = eye(homds.nphase);
        sortedevls = [];
    else
        [evc0,evl0] = eig(A0);
        [K,i] = min(abs(diag(evl0)));
        evc = evc0(:,[1:i-1 i+1:end]);%eigenvectoren zonder eigenvector horende bij 0
        evl = evl0([1:i-1 i+1:end],[1:i-1 i+1:end]);%eigenwaarden zonder 0
    
        val = 0;
        pos = 0;
        k = 1;
        for i=1:size(A0,1)-1
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
        VU = zeros(size(A0,1),size(val,2));
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
        [Q0,RU] = qr(B);
    end
    
elseif flag == 1    
    if dim == 0
        Q0 = eye(homds.nphase);
        sortedevls = [];
    else
        [evc0,evl0] = eig(A0);
        [K,i] = min(abs(diag(evl0)));
        evc = evc0(:,[1:i-1 i+1:end]);%eigenvectoren zonder eigenvector horende bij 0
        evl = evl0([1:i-1 i+1:end],[1:i-1 i+1:end]);%eigenwaarden zonder 0
    
        val = 0;
        pos = 0;
        k = 1;
        for i=1:size(A0,1)-1
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
        VU = zeros(size(A0,1),size(val,2));
        sortedevls = zeros(size(pos,2),1);
        for l = 1:size(val,2)
            VU(:,l) = evcs(:,b(l));
            sortedevls(l) = evls(b(l)); %de eigenwaarden gesorteerd van klein naar groot (enkel negatief reeel deel)
        end   
    
        VU_1 = zeros(size(A0,1),size(val,2));
        sortedevls_1 = zeros(size(pos,2),1);
        for f = 1:size(val,2)
            VU_1(:,f) = VU(:,end-f+1);
            sortedevls_1(f) = sortedevls(end-f+1);
        end
        sortedevls = sortedevls_1';
    
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
        [Q0,RU] = qr(B);
    end
end