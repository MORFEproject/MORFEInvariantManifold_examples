N=20;L=0.06;
handles=feval(@bruss);
[t,x0,options]=feval(handles{1},N);
[x1,v1]=init_EP_EP(@bruss,x0,[N;L],[2]);
opt=contset;
opt=contset(opt,'MinStepsize',1e-5);
opt=contset(opt,'MaxCorrIters',10);
opt=contset(opt,'MaxNewtonIters',20);
opt=contset(opt,'FunTolerance',1e-3);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',500);
opt=contset(opt,'Locators',[]);
[x,v,s,h]=cont(@pde_1,x1,v1,opt);
cpl(x,v,s,[41;20]);

