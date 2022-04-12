% Computes a new mesh, when given a limit cycle and new numbers of mesh
% intervals and collocation points. It will try to design the mesh such
% that the error in interpolation is spread out evenly over the whole
% cycle.
%
% =========================================================================

%function nmesh = adaptmesh(lc,ups,newntst)
function tmnew = newmesh(ups,msh,oldntst,oldncol,newntst,newncol)

% PART 1: COMPUTE ERROR IN PREVIOUS MESH-PARTITIONING
% ---------------------------------------------------

% First, 'wh' will contain the coefficients of (x-1)^oldncol
wh = [1 zeros(1,oldncol)];
for i=2:oldncol+1
    wh(2:i) = wh(1:i-1) - wh(2:i);
    wh(1) = -wh(1);
end
% Then, 'wh' will contain the coefficients of (oldncol*(x-1))^oldncol
wh = (oldncol^oldncol)*wh;

% Get the mesh interval widths 
dt = msh(2:end)-msh(1:end-1);

% The elements of sc are (1/dt)^oldncol
sc = (1./dt).^oldncol;

% Elements of hd:   
%   (the coefficients of (oldncol*(x-1)/dt)^oldncol) * (point on cycle)
% The computations are done per mesh interval. This is an approximation for
% the error made currently in this interval.
for j=1:oldntst
    hd(:,j) = sc(j)*(ups(:,(j-1)*oldncol+(1:(oldncol+1)))*wh');
end

if all(abs(hd) < 1e-7) 
    % If all errors in all intervals are close enough to zero, give default
    % value to eqf.
    cumulerror = 0:oldntst;
else
    % Make hd and dt cyclic
    hd(:,end+1) = hd(:,1);
    dt(end+1) = dt(1);

    cumulerror(1) = 0;
    pwr = 1/(oldncol+1);
    for j=1:oldntst
        % Average dt over 2 consecutive intervals (and last + first)
        averagedt = (dt(j) + dt(j+1)) / 2;
        
        % Average the error per interval hd over 2 consecutive
        % intervals and over the length (so that it is independent of
        % interval length)
        hd(:,j) = (hd(:,j+1) - hd(:,j)) / averagedt;
        
        % meshinterror = this is the best approximation of the error made
        % in a certain mesh interval (per length-unit)
        meshinterror = sum(abs(hd(:,j)).^pwr);
        
        % cumulerror = array containing cumulative errors over the mesh
        % intervals.
        cumulerror(j+1) = cumulerror(j) + dt(j)*meshinterror;
    end
end

% PART 2: FIND MESH-PARTITIONING WITH NEW NUMBER INTERVALS TO EQUALLY
% DISTRIBUTE THE ERROR DETECTED PREVIOUSLY
% -------------------------------------------------------------------

% The average error per interval with the new number of mesh intervals
newaverageerror = cumulerror(oldntst+1)/newntst;
% The optimal cumulative error-vector over the intervals (so per point)
idealcumulerror = (0:newntst)*newaverageerror;

% idealcumulerror contains the ideal cumulative error to have at each mesh
% point (over each mesh interval). We go over each of these points.
nmesh(length(idealcumulerror)) = 0;
for j=1:length(idealcumulerror)
    % At each point, find the index of the cumulative error-vector of the
    % old situation, where you just passed the ideal cumulative error...
    m = min(find(cumulerror > idealcumulerror(j)));
    if isempty(m)
        m = length(cumulerror);
    end
    % ...and go subtract 1 to get the index right before reaching the ideal
    % cumulative error
    m = m-1;
    
    % Now we will find the new mesh point by interpolating between two old
    % mesh points: the one right before, and the one right after reaching
    % the ideal cumulative error.
    t = (idealcumulerror(j) - cumulerror(m)) / (cumulerror(m+1) - cumulerror(m));
    nmesh(j) = (1-t)*msh(m) + t*msh(m+1);
end
tmnew = nmesh;
