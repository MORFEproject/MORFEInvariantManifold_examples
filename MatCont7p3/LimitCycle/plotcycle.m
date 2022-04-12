
function res = plotcycle(x,v,s,e)
%
% res = plotcycle(x,v,s,e)
%
% This function plots the cycles x. e is an array which defines 
% which coordinates of the differential equation should be used.
% This must be two or three dimensional.
%
global lds

watchon;

if length(e)==2
    p=plot(x(e(1),1),x(e(2),1),'-',x(e(1),1),x(e(2),1),'*');
    hold on;
    xmin = min(min(x((0:lds.tps-1)*lds.nphase+e(1),:)));
    xmax = max(max(x((0:lds.tps-1)*lds.nphase+e(1),:)));
    ymin = min(min(x((0:lds.tps-1)*lds.nphase+e(2),:)));
    ymax = max(max(x((0:lds.tps-1)*lds.nphase+e(2),:)));
    axis([xmin xmax ymin ymax]);
 %   xlabel('x(1)');
 %   ylabel('x(2)');
 %   pause(0.1);
    for j=1:size(x,2)
        if isempty(find([s.index]==j))
            plot(x((0:lds.tps-1)*lds.nphase+e(1),j),x((0:lds.tps-1)*lds.nphase+e(2),j),'-');
        else
            plot(x((0:lds.tps-1)*lds.nphase+e(1),j),x((0:lds.tps-1)*lds.nphase+e(2),j),'r-');
        end
    end
  %  for j=1:size(x,1)/lds.nphase
  %      plot(x((j-1)*lds.nphase+1,:),x((j-1)*lds.nphase+2,:),':');
  %  end
    hold off;
else
    p=plot3(x(e(1),1),x(e(2),1),x(e(3),1),'-');
    hold on;
    xmin = min(x(e(1),:));
    xmax = max(x(e(1),:));
    ymin = min(min(x((0:lds.tps-1)*lds.nphase+e(2),:)));
    ymax = max(max(x((0:lds.tps-1)*lds.nphase+e(2),:)));
    zmin = min(min(x((0:lds.tps-1)*lds.nphase+e(3),:)));
    zmax = max(max(x((0:lds.tps-1)*lds.nphase+e(3),:)));
    axis([xmin xmax ymin ymax zmin zmax]);
    for j=1:size(x,2)
        if isempty(find([s.index]==j))
            plot3(x(e(1),j)*ones(lds.tps,1),x((0:lds.tps-1)*lds.nphase+e(2),j),x((0:lds.tps-1)*lds.nphase+e(3),j),'-');
        else
            plot3(x(e(1),j)*ones(lds.tps,1),x((0:lds.tps-1)*lds.nphase+e(2),j),x((0:lds.tps-1)*lds.nphase+e(3),j),'r.-');
        end
    end
 %   for j=1:size(x,1)/lds.nphase
 %       plot3(x(e(1),:),x((j-1)*lds.nphase+e(2),:),x((j-1)*lds.nphase+e(3),:),':');
 %   end
    hold off;
end
res = p;
watchoff;
