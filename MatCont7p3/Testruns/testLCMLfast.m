testEquilMLfast
x1=x(1:2,s(2).index);p=[x(end,s(2).index);0.1];
[x0,v0]=init_H_LC(@MLfast,x1,p,ap1,0.0001,30,4);
opt=contset;
opt=contset(opt,'IgnoreSingularity',1);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',50);
[x2,v2,s2,h2,f2]=cont(@limitcycle,x0,v0,opt);

plotcycle(x2,v2,s2,[1 2]);
