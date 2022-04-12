function [start_point, data, msg] = GUIConnect_getPointHet(settings)
start_point = [];
msg = ''; 
data = struct();

odesys = settings.system;
x0 = settings.coord;
x1 = settings.coord_target;
x1 = x1(:);
dim = length(x0);
par = settings.parameters;

UParam1 = settings.con_UParam1;
UParam2 = settings.con_UParam2;
eps0 = settings.con_eps0;

increment = settings.option_increment;

func_handles = feval(odesys.handle);


saddle_point = x0(:);
target_saddle_point = x1(:);
data.saddle_point = saddle_point;
data.target_saddle_point = target_saddle_point;


% % % END INIT % % %
if UParam1==0 && UParam2==0
    msg = ('UParam1 or UParam2 has to be different from zero');
    return
end

if eps0 == 0
    msg = ('eps0 has to be different from zero');
    return
end

if eps0 < 0
    msg = ('eps0 has to be positive');
    return
end

check = 0;
if ~isempty(func_handles{3})
    jaco = func_handles{3};
    check = 1;
end
p = num2cell(par);

% is x_0 an equilibrium?
func = feval(func_handles{2}, 0, saddle_point, p{:});
if ~(norm(func)<1e-1)
    msg = ('start is not an equilibrium');
    return
end
% is x_1 an equilibrium?
func = feval(func_handles{2}, 0, target_saddle_point, p{:});
if ~(norm(func)<1e-1)
    msg = ('target is not an equilibrium');
    return
end

if check == 1
    A = feval(jaco, 0,saddle_point, p{:});
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

D = eig(A); %dimensie van D : dim maal 1
dim_npos = sum(real(D) > 0); 

[QU, eigvl] = computeBaseConnecHet(A,0); %eigvl = de gesorteerde eigenwaarden van klein naar groot, met positief reeel deel
%unstable

minimum = inf;
minimum2 = inf;
pos = inf;
pos2 = inf; %positie van 2e kleinste positieve eigenwaarde
%je zoekt de eigenwaarde met het kleinste positief reeel deel
for k = 1:dim
    %%%
    if real(D(k)) > 0
        if real(D(k)) < minimum
            pos2 = pos;
            minimum2 = minimum;
            minimum = real(D(k));            
            pos = k;
        elseif real(D(k)) <= minimum2            
            minimum2 = real(D(k));
            pos2 = k;       
        end        
    end
end

%dan ga je enkel c1 gebruiken, dus je kan c2 gelijk stellen aan 0
if ~(minimum == minimum2)
    UParam2 = 0;
end

if UParam1==0 && UParam2==0
    msg = ('UParam1 has to be different from zero'); return;
end

sumsquared = UParam1^2+UParam2^2;
UParam1 = UParam1/sqrt(sumsquared);
UParam2 = UParam2/sqrt(sumsquared);
start_point = saddle_point + eps0*(UParam1*QU(:,1)+UParam2*QU(:,2));

A = feval(jaco, 0,target_saddle_point, p{:});
data.A = A;
data.dim_npos = dim_npos;
data.label = 'HTHet';
data.loader = @GUIselectHTHet;




