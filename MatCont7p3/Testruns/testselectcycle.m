

OPTIONS = [];
hls = adaptx;
OPTIONS=odeset('RelTol',1e-8);
[t,y] = ode45(hls{2},[0 300],[0.3 0.5 -0.1],OPTIONS,1,0.8);

x1 = y(end,:);

[t,y] = ode45(hls{2},[0 10],x1,OPTIONS,1,0.8);

figure
plot(y(:,1),y(:,2))

p=[1;0.8];
ap=[2];

tolerance=1e-2;
[x0,v0]=initOrbLC(@adaptx,t,y,p,ap,20,4,tolerance);
opt=contset;
opt=contset(opt,'MaxNumPoints',50);
%opt=contset(opt,'TSearchOrder',0);
%opt=contset(opt,'Backward',1);
[xlcc,vlcc,slcc,hlcc,flcc]=cont(@limitcycle,x0,v0,opt);

figure
axes
plotcycle(xlcc,vlcc,slcc,[245 1 2]);