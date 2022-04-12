% Ricatti evaluation file
% =======================

function result = RicattiEvalHTHet(x0,p,unstable_flag,Y)
global HTHetds

A = cjac(HTHetds.func,HTHetds.Jacobian,x0,p,HTHetds.ActiveParams);
result = [];

if unstable_flag
    if  ~isempty(Y)
    Q0U = HTHetds.oldUnstableQ;
    % Riccati blocks from unstable eigenspace
    [U11, U12, UE21, U22] = RicattiCoeffHTHet(Q0U,A,HTHetds.npos);
    tmp = (U22*Y - Y*U11 + UE21 - Y*U12*Y)';
    for i=1:HTHetds.nphase-HTHetds.npos
        result(end+1:end+HTHetds.npos,1) = tmp(:,i);
    end
    end
else
    if  ~isempty(Y)
    Q1S = HTHetds.oldStableQ;
    % Riccati blocks from stable eigenspace
    [S11, S12, SE21, S22] = RicattiCoeffHTHet(Q1S,A,HTHetds.nneg);
    tmp = (S22*Y - Y*S11 + SE21 - Y*S12*Y)';
    for i=1:HTHetds.nphase-HTHetds.nneg
        result(end+1:end+HTHetds.nneg,1) = tmp(:,i);
    end
    end
end