function plotcycle1(x,v,s,e)
%
% plotcycle(x,v,s,e)
%
% This function plots the cycles x. e is an array whos first two
% elements define which coordinates of the differential equation 
% should be used. This must be two or three dimensional.
%

global lds

if length(e)==2
    line(x((0:lds.tps-1)*3+e(1),1),x((0:lds.tps-1)*3+e(2),1),'color','b','linestyle','none','marker','*','markerfacecolor','g','markeredgecolor','g');   
    line(x((0:lds.tps-1)*3+e(1),:),x((0:lds.tps-1)*3+e(2),:),'color','b','linestyle','-');
    line(x((0:lds.tps-1)*3+e(1),cat(1,s.index)),x((0:lds.tps-1)*3+e(2),cat(1,s.index)),'color','r','linestyle','-.','linewidth',2);
    line(x((0:lds.tps-1)*3+e(1),:)',x((0:lds.tps-1)*3+e(2),:)','color','b','linestyle',':');
    xlabel(['x(' num2str(e(1)) ')']);
    ylabel(['x(' num2str(e(2)) ')']);
else
    line(x(end*ones(lds.tps,1),:),x((0:lds.tps-1)*lds.nphase+e(1),:),x((0:lds.tps-1)*lds.nphase+e(2),:),'color','b','linestyle','-');
    line(x(end*ones(lds.tps,1),cat(1,s.index)),x((0:lds.tps-1)*lds.nphase+e(1),cat(1,s.index)),x((0:lds.tps-1)*lds.nphase+e(2),cat(1,s.index)),'color','r','linestyle','-.','linewidth',2);
    line(x(end*ones(lds.tps,1),:)',x((0:lds.tps-1)*lds.nphase+e(1),:)',x((0:lds.tps-1)*lds.nphase+e(2),:)','color','b','linestyle',':');
    view(3);
    xlabel(['par(' num2str(e(3)) ')']);
    ylabel(['x(' num2str(e(1)) ')']);
    zlabel(['x(' num2str(e(2)) ')']);
end

