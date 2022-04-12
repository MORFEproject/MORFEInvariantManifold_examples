function [tmnew,upsnew,mshnew] = newmeshcycle(ups,msh,oldntst,oldncol,newntst,newncol)

wh = 1:(oldncol+1);
wh(1) = 1;
for i=1:oldncol
  wh(i+1) = 0;
  wh((i+1):-1:2) = wh(i:-1:1) - wh((i+1):-1:2);
  wh(1) = -wh(1);
end
wh = (oldncol^oldncol)*wh;
dt = msh(2:(oldntst+1))-msh(1:oldntst);
k=find(dt==0);
if ~isempty(k)
    msh(k)=[];
    ups(:,k)=[];
    oldntst = oldntst-size(k,2);
    dt = msh(2:(oldntst+1))-msh(1:oldntst);
end
sc = (1./dt).^oldncol;
for j=1:oldntst
  hd(:,j) = sc(j)*(ups(:,(j-1)*oldncol+(1:(oldncol+1)))*wh');
end
if all(abs(hd)<1e-7),eqf = 0:oldntst;
else
  hd(:,oldntst+1) = hd(:,1);

  pwr = 1/(oldncol+1);
  eqf(1) = 0;
  for j=1:oldntst
    if j<oldntst
      dtav = (dt(j+1)+dt(j))/2;
    else
      dtav = (dt(1)+dt(j))/2;
    end
    hd(:,j) = (1/dtav)*(hd(:,j+1)-hd(:,j));
    e = sum(abs(hd(:,j)).^pwr);
    eqf(j+1)=eqf(j)+dt(j)*e;

  end
end;

% make new time mesh such that eqf is more or less equally
% distributed over the interval
dal = eqf(oldntst+1)/newntst;
uneq = (0:newntst)*dal;

ial = ordr(eqf,uneq);
tmnew = zeros(1,newntst+1);
for j=1:(newntst+1)
  k = ial(j);
  t = (uneq(j)-eqf(k))/(eqf(k+1)-eqf(k));
  tmnew(j) = (1-t)*msh(k)+t*msh(k+1);
end
mshnew = msh;
upsnew  = ups;