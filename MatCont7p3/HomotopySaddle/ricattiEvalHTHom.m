% Ricatti evaluation file
% =======================

function result = ricattiEvalHTHom(x0,p,YU,YS)
global HTHomds

Q0S = HTHomds.oldStableQ;
Q0U = HTHomds.oldUnstableQ;
A = cjac(HTHomds.func,HTHomds.Jacobian,x0,p,HTHomds.ActiveParams);
result = [];

if ~isempty(YU)
    % Riccati blocks from unstable eigenspace
    [U11, U12, UE21, U22] = ricattiCoeff(Q0U,A,HTHomds.npos);
    tmp = (U22*YU - YU*U11 + UE21 - YU*U12*YU)';
    for i=1:HTHomds.nneg
        result(end+1:end+HTHomds.npos,1) = tmp(:,i);
    end
end

if ~isempty(YS)
    % Riccati blocks from stable eigenspace
    [S11, S12, SE21, S22] = ricattiCoeff(Q0S,A,HTHomds.nneg);
    tmp = (S22*YS - YS*S11 + SE21 - YS*S12*YS)';
    for i=1:HTHomds.npos
        result(end+1:end+HTHomds.nneg,1) = tmp(:,i);
    end
end
