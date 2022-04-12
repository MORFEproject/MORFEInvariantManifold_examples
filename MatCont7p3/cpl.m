function varargout = cplanim(x,v,s,varargin)
%
% plots computed curve: cpl(x,v,s)
% dim(x) = 2 or 3
%
if nargin > 3
  x = x(varargin{1},:);
  v = v(varargin{1},:);
end
[dim,stop]=size(x);

watchon;
switch dim
  case 2
    res=cpl2d(x,v,s);
  case 3
    res=cpl3d(x,v,s);
  otherwise
    error('This dimension is not supported');
end
watchoff;
function res = cpl2d(x,v,s)
line(x(1,:),x(2,:),'linestyle','-','color','b');

xindex=cat(1,s.index);
xindex=xindex(2:end-1);
line(x(1,xindex),x(2,xindex),'linestyle','none', 'marker', '.','color','r');
%%line(x(1,cat(1,s.index)),x(2,cat(1,s.index)),'linestyle','.','color','r'); %%
d=axis;
if size(s,1)~=2
    skew = 0.03*[d(2)-d(1) d(4)-d(3)];
%    xindex=cat(1,s.index);
%    xindex=xindex(2:end-1);
    s(1).label=''; s(end).label=''; 
    text( x(1,xindex)+skew(1),x(2,xindex)+skew(2), {s(2:end-1).label});
end
res =2;

function res = cpl3d(x,v,s)

view(3);
line(x(1,:),x(2,:),x(3,:),'linestyle','-','color','b');
xindex=cat(1,s.index);
xindex=xindex(2:end-1);
line(x(1,xindex),x(2,xindex),x(3,xindex),'linestyle','none', 'marker', '.','color','r');
%line(x(1,cat(1,s.index)),x(2,cat(1,s.index)),x(3,cat(1,s.index)),'linestyle','none', 'marker', '.','color','r'); %*
d=axis;
if size(s,1)~=2
    skew = 0.03*[d(2)-d(1) d(4)-d(3) d(6)-d(5)];
 %   xindex=cat(1,s.index);
 %   xindex=xindex(2:end-1);
    s(1).label=''; s(end).label='';
    text( x(1,xindex)+skew(1),x(2,xindex)+skew(2), x(3,xindex)+skew(3) ,{s(2:end-1).label});
   
end
res =3;



%SD:plots the 2d or 3d curve
