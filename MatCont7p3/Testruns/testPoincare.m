
TSTP=@testEV;
OPTIONS = odeset('RelTol',1e-8,'Events',TSTP);
%OPTIONS=[]
hls = adaptx;
[t,y,TE,YE,IE] = ode45(hls{2},[0 300],[0.3 0.5 -0.1],OPTIONS,1,0.8);
x0 = y(end,:);
[t,y,TE,YE,IE] = ode45(hls{2},[0 10],x0,OPTIONS,1,0.8);

figure
plot(y(:,1),y(:,2))
TE,YE,IE

function [value,isterminal,direction]= testEV(t,y,varargin)
value=[y(1)-0.2;y(2)-0.3];
isterminal=zeros(2,1);
direction=ones(2,1);
end