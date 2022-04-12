% Ricatti evaluation file
% =======================

function result = HSN_ricattiEval(x0,p,YU,YS)

global homds

Q0S = homds.oldStableQ;
Q0U = homds.oldUnstableQ;
A = cjac(homds.func,homds.Jacobian,x0,p,homds.ActiveParams);
result = [];

if ~isempty(YU)
    % Riccati blocks from unstable eigenspace
    [U11, U12, UE21, U22] = ricattiCoeff(Q0U,A,homds.npos);
    tmp = (U22*YU - YU*U11 + UE21 - YU*U12*YU)';
    for i=1:homds.nneg+1
        result(end+1:end+homds.npos,1) = tmp(:,i);
    end
end

if ~isempty(YS)
    % Riccati blocks from stable eigenspace
    [S11, S12, SE21, S22] = ricattiCoeff(Q0S,A,homds.nneg);
    tmp = (S22*YS - YS*S11 + SE21 - YS*S12*YS)';
    for i=1:homds.npos+1
        result(end+1:end+homds.nneg,1) = tmp(:,i);
    end
end

