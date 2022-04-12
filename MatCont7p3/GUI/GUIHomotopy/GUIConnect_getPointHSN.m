function [start_point, data, msg] = GUIConnect_getPointHSN(settings)
start_point = [];
msg = ''; 
data = struct();

odesys = settings.system;
x0 = settings.coord;
dim = length(x0);
par = settings.parameters;

UParam1 = settings.con_UParam1;
eps0 = settings.con_eps0;

increment = settings.option_increment;

func_handles = feval(odesys.handle);


saddle_point = x0(:);
data.saddle_point = saddle_point;

% % % END INIT % % %


if UParam1==0
    msg = 'UParam1 has to be different from zero';
    return;
end

if eps0 == 0
    msg = 'eps0 has to be different from zero';
    return;
end

if eps0 < 0
    msg = 'eps0 has to be positive';
    return;
end

check = 0;
if ~isempty(func_handles{3})
    jaco = func_handles{3};
    check = 1;
end
p = num2cell(par);

if check ==1
    A = feval(jaco, 0, saddle_point, p{:});
else
    point = saddle_point;
    x1 = point;
    x2 = point;
    for i=1:dim
        x1(i) = x1(i)-increment;
        x2(i) = x2(i)+increment;
        j(:,i) = feval(func_handles{2}, 0, x2, p{:})-feval(func_handles{2}, 0, x1, p{:});
        x1(i) = point(i);
        x2(i) = point(i);
    end
    A = j/(2*increment);
end

[V,D] = eig(A); %dimensie van D : dim maal 1
[K,i] = min(abs(diag(D)));%eigenwaarde 0
eig0 = V(:,i); 

% is x_0 an equilibrium?
func = feval(func_handles{2}, 0, saddle_point, p{:});

if ~(norm(func)<=1e-1)
    msg = ('The starting point is not an equilibrium');   
    return
end

if ~(abs(K)<=1e-1)
    msg = ('The starting point is not a saddle-node');   
    return
end    

B = D([1:i-1 i+1:end],[1:i-1 i+1:end]);%je sluit eigenwaarde 0 uit
dim_nneg = sum(real(diag(B)) < 0);
if (dim_nneg == dim-1)   
    if min(abs(real(diag(B)))) < 1e-10
        dim_nneg = dim_nneg -1;    
    end
end
if (dim_nneg == 0)
    if min(abs(real(diag(B)))) < 1e-10
        dim_nneg = dim_nneg +1;
    end
end
dim_npos = dim-dim_nneg-1; 

UParam1 = UParam1/abs(UParam1);
start_point = saddle_point + eps0*UParam1*eig0;


data.A = A;
data.dim_npos = dim_npos;
data.label = 'HTHSN';
data.loader = @GUIselectHTHSN;
data.eig0 = eig0;

