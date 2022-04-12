[x0,v0]=init_EP_EP(@adaptx,[0;0;0],[-10;1],[1]);
opt=contset;opt=contset(opt,'Singularities',1);
[x,v,s,h,f]=cont(@equilibrium,x0,[],opt);

x1=x(1:3,s(2).index);p=[x(end,s(2).index);1];
[x0,v0]=init_H_LC(@adaptx,x1,p,[1],1e-6,20,4);
opt = contset(opt,'MaxNumPoints',10);
opt = contset(opt,'Multipliers',1);
opt = contset(opt,'Adapt',1);
[xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt);
par=slc(end).data.parametervalues;
[x1,v1] = init_LC_LC(@adaptx, xlc, vlc, slc(end), par, 1, 20, 4);
opt = contset(opt,'PRC',1);
opt = contset(opt,'dPRC',1);
opt = contset(opt,'Input',1);
opt = contset(opt,'MaxNumPoints',20);
[xlc1,vlc1,slc1,hlc1,flc1]=cont(@limitcycle,x1,v1,opt);

fvector=flc1(:,20);
plot(fvector(22:102),'r')
hold on
plot(fvector(103:183),'b')
