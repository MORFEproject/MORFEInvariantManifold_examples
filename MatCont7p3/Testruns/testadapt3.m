function testadapt3
global xlc vlc slc opt xpd vpd spd hpd fpd
testadapt1;
disp('>> [x0,v0]=init_PD_PD(@adaptx,xlc,slc(4),[1 2],20,4);');
[x0,v0]=init_PD_PD(@adaptx,xlc,slc(4),[1 2],20,4);
disp('>> opt = contset; opt = contset(opt,''Singularities'',1);');
opt=contset;opt=contset(opt,'Singularities',1);
disp('>> [xpd,vpd,spd,hpd,fpd]=cont(@perioddoubling,x0,v0,opt);');
[xpd,vpd,spd,hpd,fpd]=cont(@perioddoubling,x0,v0,opt);
figure;
axes;
disp('>> cpl(xpd,vpd,spd,[245 246]);');
cpl(xpd,vpd,spd,[245 246]);
