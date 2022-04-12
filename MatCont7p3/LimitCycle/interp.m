% Re-interpolates a given vector (f), defined on oldmsh, on the new
% mesh
%
% =========================================================================

function newvect = interp(oldmsh,oldncol,f,newmsh,newncol)

% Length of a small subinterval: between 2 finemesh points
subintdt = (newmsh(2:end)-newmsh(1:(end-1)))/newncol;
newlen = length(newmsh);
oldlen = length(oldmsh);
oldsubintdt = (oldmsh(2:end) - oldmsh(1:end-1)) / oldncol;

for i=1:newncol
    % Tempmesh is a mesh containing the positions of the collocation point
    % discussed now, in each subinterval
    tempmesh = newmsh(1:(newlen-1)) + (i-1)*subintdt;
    
    for j=1:(newlen-1)
        % Get the index of the old meshpoint that is the last one before
        % reaching the position of the new fine meshpoint
        m = min(find(oldmsh > tempmesh(j)));
        if isempty(m)
            m = length(oldmsh);
        end
        m = m-1;
        
        % First, x contains the positions of the old fine meshpoints within
        % the mesh interval in which the fine meshpoint will be that is
        % studied right now.
        x = oldmsh(m) + (0:oldncol)*oldsubintdt(m);
        % Then x contains the difference in position between x's previous
        % elements, and the current fine meshpoint.
        x = tempmesh(j) - x;
        
        % Store the found differences in a matrix, each column in a
        % different order.
        xmat = x(reshape(rem((((oldncol+1)^2)-1:-1:oldncol+1),oldncol+1)+1,...
                        oldncol,oldncol+1));
        % Compute the weights for interpolation
        w = prod(xmat,1)./prod(xmat-x(ones(1,oldncol),:),1);
    
        % Assign new value to fine meshpoint.
        % The new value is found by interpolating all fine meshpoints from
        % the old meshconfiguration, using w as weights.
        newvect(:,(j-1)*newncol+i) = f(:,(m-1)*oldncol+(1:(oldncol+1)))*w';
    end
end
% Assign value to the last collocation point, the last point of the cycle.
% Logically, its coordinates remain the same as before.
newvect(:,(newlen-1)*newncol+1) = f(:,(oldlen-1)*oldncol+1);