testEquilMLfast
x1=x(1:2,s(2).index);p=[x(end,s(2).index);0.1];
[x0,v0]=init_H_LC(@MLfast,x1,p,ap1,0.0001,30,4);
opt=contset;
opt=contset(opt,'IgnoreSingularity',1);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',50);
opt=contset(opt,'FunTolerance',0.0000001);
opt=contset(opt,'VarTolerance',0.0000001);
[x2,v2,s2,h2,f2]=cont(@limitcycle,x0,v0,opt);
[x0,v0]=init_LPC_LPC(@MLfast,x2,s2(2),[1 2],30,4);
opt=contset;
opt=contset(opt,'FunTolerance',0.0001);
opt=contset(opt,'VarTolerance',0.0001);
opt=contset(opt,'MaxNumPoints',30);
%opt=contset(opt,'Backward',1);
opt=contset(opt,'Singularities',1);
[x3,v3,s3,h3,f3]=cont(@limitpointcycle,x0,v0,opt);
plotcycle(x3,v3,s3,[1 2]);