%
% PRC = calcPRC(x0, Input, which)
% Computes the phase response curve of a limit cycle and/or its derivative.
% x0 = current limit cycle with attached parameters in 1 vector
% Input = pulse size (vector or double value)
% which = [1 0] for PRC, [0 1] for dPRC or [1 1] for both

function [PRC, dPRC] = calcPRC(x0, Input, which)
global lds cds

dimX = lds.tps * lds.nphase;
x = x0(1:dimX);
tmpcycle = reshape(x,lds.nphase,lds.tps);
if length(lds.ActiveParams) == 1
    newperiod = x0(dimX+1);
    newpars = x0(dimX+2);
else
    newperiod = lds.T;
    newpars = x0(dimX+1:dimX+length(lds.ActiveParams));
end
p = lds.P0;
p(lds.ActiveParams) = newpars;

pos = (gl_pos(lds.ncol) + 1)/2;
wcol = gl_weight(lds.ncol)/2;

% Compose reduced big jacobian
% ----------------------------
% Remove derivative by parameter
if (~isfield(cds,'adapted') || cds.adapted == 0) && ~isempty(cds.oldJac) && size(cds.oldJac,1) == (lds.tps*lds.nphase +1)
    RJ = cds.oldJac;
else
    RJ = cjac(cds.curve_func,cds.curve_jacobian,x0(:,end),[]);
end
RJ = RJ(:,1:end-1);
% The resulting matrix is square and non-singular
% 
% size(RJ)
% size([zeros(1,lds.tps*lds.nphase) -1/newperiod])
vltemp = [zeros(1,lds.tps*lds.nphase) -1/newperiod] / RJ;
vl = vltemp(1:end-1)';

if isa(Input, 'function_handle')
    [val,ind] = max(tmpcycle(1,:));
    
    % Remove vl(0)
    vl(end-lds.nphase+1:end,:) = [];
    
    colmesh = [];    
    range1 = 1:(lds.ncol + 1);
    for i=1:lds.ntst    
        % Collocation mesh
        colmesh(range1(1:end-1)) = lds.msh(i) + lds.dt(i) * pos;
        range1 = range1 + lds.ncol;
    end
    inpmesh = colmesh;
    inpmesh(end+1:end+ind-1) = inpmesh(1:ind-1) + 1;
    inpmesh = inpmesh(ind:end) - inpmesh(ind);
        
    tmpcycle(:,end+1:end+ind-1) = tmpcycle(:,1:ind-1);
    tmpcycle = tmpcycle(:,ind:end);
    
    % Reshape vl, so that each column is the vector in one point
    vl = reshape(vl,lds.nphase,lds.tps-1);
    vl(:,end+1:end+ind-1) = vl(:,1:ind-1);
    vl = vl(:,ind:end);
    
    for i=1:length(vl(1,:))
        dvl(:,i) = - newperiod * cjac(lds.func,lds.Jacobian,tmpcycle(:,i),p,lds.ActiveParams)' * vl(:,i);
    end
    
    range1 = 1:(lds.ncol + 1);
    for i=1:lds.ntst      
        % collocation point coordinates
        xval = tmpcycle(:,range1) * lds.wt;
        
        border = lds.ncol;
        if i == lds.ntst
            border = border +1;
        end
        
        for j=1:border
            % Evaluate function, pass limit cycle, coordinates and mesh value
            Inputv = feval(Input,xval(:,j));
            if size(Inputv,2) > 1
                Input = Inputv';
            end
            PRC(range1-1+j) = vl(:,range1-1+j)' * Inputv;    
            dPRC(rang1-1+j) = dvl(:,range1-1+j)' * Inputv;
        end 

        range1 = range1 + lds.ncol;
    end    
    dPRC = dPRC/newperiod;
    
    % Plot PRC
    if which(1) ~= 0
        figure(which(1));
        plot(inpmesh, PRC,'b','LineWidth',1);
        hold on
        plot(inpmesh,zeros(size(inpmesh)),'k','LineWidth',1);
        %%plot voltage on top:
        %plot(inpmesh*newperiod,tmpcycle(1,:)*50,'r:','LineWidth',1);
        xlabel('Phase')
        ylabel('Response')
        title('PRC')
    end  
    % Plot dPRC
    if which(2) ~= 0
        figure(which(2));     
        plot(inpmesh, dPRC,'r','LineWidth',1);    
        hold on
        plot(inpmesh,zeros(size(inpmesh)),'k','LineWidth',1);
        xlabel('Phase')
        ylabel('Derivative of response')
        title('dPRC')
    end
    
elseif isa(Input, 'double')    

    if length(Input) == 1
        Input = [Input; zeros(lds.nphase-1,1)];
    end
    
    if length(Input) ~= lds.nphase
        error('Input requires a vector with length = number of state variables');
    end
    
    % Reshape vl, so that each column is the vector in one point
    vl = reshape(vl,lds.nphase,lds.tps);
    
    % Rescale v-values:
    range1 = 1:lds.ncol;
    for i=1:lds.ntst  
        % Divide v-values in collocation points by Gauss-Legendre weights
        % and mesh-interval lengths
        vl(:,range1) = ( vl(:,range1) ./ repmat(wcol',lds.nphase,1) ) / lds.dt(i);
        range1 = range1 + lds.ncol;
    end
    
    % Transform values in collocation points to values in mesh points
    % mesh2colloc is a block-structured matrix
    mesh2colloc = sparse(lds.tps,lds.tps);
    for i=1:((lds.tps-1)/lds.ncol)
        mesh2colloc((i-1)*lds.ncol+1:i*lds.ncol+1,(i-1)*lds.ncol+1:i*lds.ncol) = lds.wt;
    end
    mesh2colloc(1,end) = 1/2;
    mesh2colloc(end,end) = 1/2;

    vl = vl / mesh2colloc;
    p = num2cell(p);
    
    for i=1:length(vl(1,:))
        dvl(:,i) = - newperiod * cjac(lds.func,lds.Jacobian,tmpcycle(:,i),p,lds.ActiveParams)' * vl(:,i);
    end
   
    % Compute Inputv
    if size(Input,2) > 1
        Input = Input';
    end
    
    % Compute PRC
    range1 = 1:(lds.ncol + 1);
    for i=1:lds.ntst    
        border = lds.ncol;
        if i == lds.ntst
            border = border +1;
        end
        % For whole cycle (all mesh points), compute Iv(t), based on vl(t)
        PRC(range1(1:border)) = vl(:,range1(1:border))' * Input;
        dPRC(range1(1:border)) = dvl(:,range1(1:border))' * Input;
    
        range1 = range1 + lds.ncol;
    end
    inpmesh = lds.finemsh;
        
    % Make the PRC and collocation mesh start at the spike-top
    [val,ind] = max(tmpcycle(1,:));
    PRC(end+1:end+ind-1) = PRC(1:ind-1);
    PRC = PRC(ind:end);
    dPRC(end+1:end+ind-1) = dPRC(1:ind-1);
    dPRC = dPRC(ind:end);
    dPRC = dPRC/newperiod;
    
    tmpcycle(:,end+1:end+ind-1) = tmpcycle(:,1:ind-1);
    tmpcycle = tmpcycle(:,ind:end);

    inpmesh(end+1:end+ind-1) = inpmesh(1:ind-1) + 1;
    inpmesh = inpmesh(ind:end) - inpmesh(ind);    

    % Plot PRC
    if which(1) ~= 0
%         figure(MC.PRC(end));
        figure(which(1));
%         plot(inpmesh*newperiod, PRC,'r','LineWidth',1);
        plot(inpmesh, PRC,'b','LineWidth',1);
        hold on
        plot(inpmesh,zeros(size(inpmesh)),'k','LineWidth',1);
        %%plot voltage on top:
        %plot(inpmesh*newperiod,tmpcycle(1,:)*50,'r:','LineWidth',1);
        xlabel('Phase')
        ylabel('Response')
        title('PRC')
    end  
    % Plot dPRC
    if which(2) ~= 0
%         figure(MC.dPRC(end));
        figure(which(2));    
%         plot(inpmesh*newperiod, dPRC,'g','LineWidth',1);   
        plot(inpmesh, dPRC,'r','LineWidth',1);    
        hold on
        plot(inpmesh,zeros(size(inpmesh)),'k','LineWidth',1);
        xlabel('Phase')
        ylabel('Derivative of response')
        title('dPRC')
    end
    
else
    error('Wrong type of Input given.');
end

