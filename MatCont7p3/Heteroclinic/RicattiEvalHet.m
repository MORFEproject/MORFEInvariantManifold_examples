% Ricatti evaluation file
% =======================

function result = RicattiEvalHet(x0,p,unstable_flag,Y)
global hetds

A = cjac(hetds.func,hetds.Jacobian,x0,p,hetds.ActiveParams);
result = [];

if unstable_flag
    if  ~isempty(Y)
    Q0U = hetds.oldUnstableQ;
    % Riccati blocks from unstable eigenspace
    [U11, U12, UE21, U22] = RicattiCoeffHet(Q0U,A,hetds.npos);
    tmp = (U22*Y - Y*U11 + UE21 - Y*U12*Y)';
    for i=1:hetds.nphase-hetds.npos
        result(end+1:end+hetds.npos,1) = tmp(:,i);
    end
    end
else
    if  ~isempty(Y)
    Q1S = hetds.oldStableQ;
    % Riccati blocks from stable eigenspace
    [S11, S12, SE21, S22] = RicattiCoeffHet(Q1S,A,hetds.nneg);
    tmp = (S22*Y - Y*S11 + SE21 - Y*S12*Y)';
    for i=1:hetds.nphase-hetds.nneg
        result(end+1:end+hetds.nneg,1) = tmp(:,i);
    end
    end
end