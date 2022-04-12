function  [Q0,sortedevls] =  computeBaseHTHSN(A0,flag,dim)
%flag == 0: unstable, flag ==1: stable
global HTHSNds

if flag == 0
    if dim == 0
        Q0 = eye(HTHSNds.nphase);
        sortedevls = [];
    else
        [evc0,evl0] = eig(A0);
        [K,i] = min(abs(diag(evl0)));
        evc = evc0(:,[1:i-1 i+1:end]);%eigenvectoren zonder eigenvector horende bij 0
        evl = evl0([1:i-1 i+1:end],[1:i-1 i+1:end]);%eigenwaarden zonder 0
            
        evl = diag(evl);
        pos = find(real(evl)>0);
        evlsr = real(evl(pos));
        evls = evl(pos);
        evcs = evc(:,pos); % the eigenvectors involved
        
        [a,b] = sort(evlsr); % result: a(i) = val(b(i))
        VU = evcs(:,b);
        sortedevls = evls(b);                 
        
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
        Q0 = eye(HTHSNds.nphase);
        sortedevls = [];
    else
        [evc0,evl0] = eig(A0);
        [K,i] = min(abs(diag(evl0)));
        evc = evc0(:,[1:i-1 i+1:end]);%eigenvectoren zonder eigenvector horende bij 0
        evl = evl0([1:i-1 i+1:end],[1:i-1 i+1:end]);%eigenwaarden zonder 0
    
        evl = diag(evl);
        pos = find(real(evl)<0);
        evlsr = real(evl(pos));
        evls = evl(pos);
        evcs = evc(:,pos); % the eigenvectors involved
        
        [a,b] = sort(evlsr); % result: a(i) = val(b(i))
        VU = evcs(:,b);
        sortedevls = evls(b);     
                  
    
        VU_1 = zeros(size(A0,1),size(VU,2));
        sortedevls_1 = zeros(size(pos,2),1);
        for f = 1:size(VU,2)
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