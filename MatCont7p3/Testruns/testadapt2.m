function testadapt2
global x v s h f xlc slc vlc xlc2 vlc2 slc2 hlc2 flc2 opt
testadapt1;clc;clf;
disp('>> [x1,v1]=init_PD_LC(@adaptx,xlc,slc(4),40,4,1e-6);');
[x1,v1]=init_PD_LC(@adaptx,xlc,slc(4),40,4,1e-6);
disp('>> opt=contset(opt,''MaxNumPoints'',250);');
opt = contset(opt,'MaxNumPoints',250);
disp('>> [xlc2,vlc2,slc2,hlc2,flc2]=cont(@limitcycle,x1,v1,opt);');
[xlc2,vlc2,slc2,hlc2,flc2]=cont(@limitcycle,x1,v1,opt);
disp('>> plotcycle(xlc2,vlc2,slc2,[size(xlc2,1) 1 2]);');
figure
axes
plotcycle(xlc2,vlc2,slc2,[size(xlc2,1) 1 2]);
