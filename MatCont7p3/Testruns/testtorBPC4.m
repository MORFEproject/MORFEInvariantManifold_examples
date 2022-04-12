p=[0.5;-0.6;0.6;0.32858;0.93358;-0.9;0];
[x0,v0]=init_EP_EP(@torBPC,[0;0;0],p,[6]);
opt=contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',50);
[x,v,s,h,f]=cont(@equilibrium,x0,[],opt);