global cds
p=[0];ap=[1];
[x0,v0]=init_EP_EP(@bratu,[0;0],p,ap);
opt=contset;
opt=contset(opt,'MaxNumPoints',50);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Userfunctions',1);
UserInfo.name='userf1';
UserInfo.state=1;
UserInfo.label='u1';
opt=contset(opt,'UserfunctionsInfo',UserInfo);
[x,v,s,h,f]=cont(@equilibrium,x0,[],opt);
[x,v,s,h,f]=cont(x,v,s,h,f,cds);
cpl(x,v,s,[3 1 2]);

