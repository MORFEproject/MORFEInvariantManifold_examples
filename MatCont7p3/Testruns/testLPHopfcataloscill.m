p=[2.5;2.204678;10;0.0675;1;0.1;0.4];
ap1=[2];
[x0,v0]=init_EP_EP(@cataloscill,[0.001137;0.891483;0.062345],p,ap1);
opt=contset;
opt=contset(opt,'MaxStepSize',0.025);
opt=contset(opt,'MaxNumPoints',78);
opt=contset(opt,'Singularities',1);
[x,v,s,h,f]=cont(@equilibrium,x0,[],opt);
cpl(x,v,s,[4,1]);
'press any key'
pause
x1=x(1:3,s(3).index);
p(ap1)=x(end,s(3).index);
[x0,v0]=init_LP_LP(@cataloscill,x1,p,[2,7]);
[x2,v2,s2,h2,f2]=cont(@limitpoint,x0,v0,opt);
hold on;
cpl(x2,v2,s2,[4,1]);
'press any key'
pause
opt=contset(opt,'Backward',1);
[x3,v3,s3,h3,f3]=cont(@limitpoint,x0,v0,opt);
hold on
cpl(x3,v3,s3,[4,1]);
'press any key'
pause
x1=x(1:3,s(2).index);
p(ap1)=x(end,s(2).index);
[x0,v0]=init_H_H(@cataloscill,x1,p,[2 7]);
opt=contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxStepsize',0.1);
opt=contset(opt,'backward',0);
[x4,v4,s4,h4,f4]=cont(@hopf,x0,v0,opt);

hold on;
cpl(x4,v4,s4,[4 1]);
