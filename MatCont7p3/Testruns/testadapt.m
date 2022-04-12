function testadapt
global x v s h f opt
disp('>> init;');

disp('>> [x0,v0]=init_EP_EP(@adaptx,[0;0;0],[-10;1],[1]);');
disp('>> opt = contset; opt = contset(opt,''Singularities'',1);');
opt=contset;opt=contset(opt,'Singularities',1);
[x0,v0]=init_EP_EP(@adaptx,[0;0;0],[-10;1],[1]);
disp('>> [x,v,s,h,f]=cont(@equilibrium,x0,[],opt);');
[x,v,s,h,f]=cont(@equilibrium,x0,[],opt);
