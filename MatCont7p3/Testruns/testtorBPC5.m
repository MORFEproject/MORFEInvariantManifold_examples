clear all
p=[0.5;-0.6;0.6;0.32858;0.93358;-0.9;0];
[x0,v0]=init_EP_EP(@torBPC,[0;0;0],p,[6]);
opt=contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',50);
[x,v,s,h,f]=cont(@equilibrium,x0,[],opt);
x1=x(1:3,s(2).index);p(6)=x(end,s(2).index);ap=6;
[x0,v0]=init_H_LC(@torBPC,x1,p,ap,0.0001,25,4);
opt=contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',150);
opt=contset(opt,'Multipliers',1);
opt=contset(opt,'Adapt',5);
[xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt);
plotcycle(xlc,vlc,slc,[1 2]);