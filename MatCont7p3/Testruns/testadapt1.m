function testadapt1
clear all
global x v s h f xlc vlc slc hlc flc opt
testadapt;
disp('>> x1=x(1:3,s(2).index);p=[x(end,s(2).index);1];')
x1=x(1:3,s(2).index);p=[x(end,s(2).index);1];
disp('>> [x0,v0]=init_H_LC(@adaptx,x1,p,[1],1e-6,20,4);');
[x0,v0]=init_H_LC(@adaptx,x1,p,[1],1e-6,20,4);
disp('>> opt = contset(opt,''MaxNumPoints'',200);');
opt = contset(opt,'MaxNumPoints',200);
disp('>> opt = contset(opt,''Multipliers'',1);');
opt = contset(opt,'Multipliers',1);
disp('>> opt = contset(opt,''Adapt'',1);');
opt = contset(opt,'Adapt',1);
disp('>> [xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt);');
[xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt);
disp('>> plotcycle(xlc,vlc,slc,[size(xlc,1) 1 2]);');
figure
axes
plotcycle(xlc,vlc,slc,[size(xlc,1) 1 2]);
% [x0,v0]=init_LC_LC(@adaptx,xlc,vlc,slc(end),[1 2],20,4);
% %disp('>> [xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt);');
% [xlcc,vlcc,slcc,hlcc,flcc]=cont(@limitcycle,x0,v0,opt);
% disp('>> plotcycle(xlcc,vlcc,slcc,[size(xlcc,1) 1 2]);');
% figure
% axes
% plotcycle(xlcc,vlcc,slcc,[size(xlcc,1) 1 2]);
