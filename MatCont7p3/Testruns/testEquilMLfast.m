clear all
p=[0.11047;0.1];ap1=[1];
[x0,v0]=init_EP_EP(@MLfast,[0.047222;0.32564],p,ap1);
opt=contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',65);
opt=contset(opt,'MinStepsize',0.00001);
opt=contset(opt,'MaxStepsize',0.01);
opt=contset(opt,'Backward',1);

[x,v,s,h,f]=cont(@equilibrium,x0,[],opt);
cpl(x,v,s,[3 1]);
