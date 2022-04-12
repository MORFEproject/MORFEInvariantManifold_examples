% Ricatti evaluation file
% =======================

function result = ricattiEvalHTHSN(x0,p,YU,YS)
global HTHSNds

Q0S = HTHSNds.oldStableQ;
Q0U = HTHSNds.oldUnstableQ;
A = cjac(HTHSNds.func,HTHSNds.Jacobian,x0,p,HTHSNds.ActiveParams);
result = [];

if ~isempty(YU)
    % Riccati blocks from unstable eigenspace
    [U11, U12, UE21, U22] = ricattiCoeff(Q0U,A,HTHSNds.npos);
    tmp = (U22*YU - YU*U11 + UE21 - YU*U12*YU)';
    for i=1:HTHSNds.nneg+1
        result(end+1:end+HTHSNds.npos,1) = tmp(:,i);%eerst eerste rij, dan 2de, enzo
    end
end

if ~isempty(YS)
    % Riccati blocks from stable eigenspace
    [S11, S12, SE21, S22] = ricattiCoeff(Q0S,A,HTHSNds.nneg);
    tmp = (S22*YS - YS*S11 + SE21 - YS*S12*YS)';
    for i=1:HTHSNds.npos+1
        result(end+1:end+HTHSNds.nneg,1) = tmp(:,i);
    end
end
