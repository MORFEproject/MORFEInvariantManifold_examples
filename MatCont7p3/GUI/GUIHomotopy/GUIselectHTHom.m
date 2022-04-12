function [newsettings, msg] = GUIselectHTHom(solution)
newsettings = []; msg = '';

settings = solution.settings;
HTHomds = struct();
HTHomds.gds = struct();
saddle_point = solution.connectData.saddle_point;
dim_npos = solution.connectData.dim_npos;
QU1 = computeBaseConnecHom(solution.connectData.A, 1);
dim = length(saddle_point);


prompt  = {'Enter the number of test intervals';'Enter the number of collocation points'};
title   = 'Choose ntst and ncol';
lines   = 1;
def     = {'40','4'};
answer  = inputdlg(prompt,title,lines,def);
if isempty(answer)||isempty(str2double(answer))
    ntst=40;
    ncol=4;
else
    ntst=str2double(answer{1});
    ncol=str2double(answer{2});    
end

y1 = [solution.y, zeros(size(solution.y, 1), 1)];
  
for k = 1:size(solution.y,1)         
    y1(k,dim+1) = norm(saddle_point'-solution.y(k,:));        
end

x = y1';

ind = size(x,2); %number of points
while ((ind > 1) && isnan(x(end,ind))) || ((ind > 1) && (x(end,ind) - x(end,ind-1) >= 0))
    ind = ind-1;
end

%ind = 0: monotonous increasing
if ind <= 10
    %we take the whole orbit as starting orbit
    ind = size(x,2);
end
if ind == size(x,2)
    %monotonous decreasing
    verschil = x(end,2:end) - x(end,1:end-1);
    pos = find(verschil > 0);
    if size(pos,2) == 0
        msg = ('No enough data available!');
        return;
    end
end

x=x(1:dim,1:ind);%omit eps1, just keep the coordinates
t=solution.t(1:ind)';
tn=(t-t(1))/(t(end)-t(1));% tn(1) = 0 and tn(end)=1



[a,x,tn] = newmeshHT(x,tn,size(x,2)-1,1,ntst,ncol);

x=interp(tn,1,x,a,ncol);

HTHomds.msh = a;
s=struct();
s(1).data.timemesh = a;
s(1).data.ntst = ntst;
s(1).data.ncol = ncol;
s(1).data.parametervalues = solution.param;
epsilon1 = norm(x(:,end)-saddle_point);

HTHomds.SParams = [];
for i = 1:dim_npos
    HTHomds.SParams{i,1} = strcat('SParam',num2str(i));    
    HTHomds.SParams{i,2} = 1/epsilon1*(x(:,end)-saddle_point)'*QU1(:,end-dim_npos+i);
end


HTHomds.P0 = solution.param;
HTHomds.ntst = ntst;
HTHomds.ncol = ncol;
HTHomds.nphase = dim;
x = reshape(x,size(x,2)*size(x,1),1);
HTHomds.T = (t(end)-t(1))/2;%half of the period
HTHomds.eps0 = settings.con_eps0;
HTHomds.eps1 = epsilon1;
s(1).data.T = HTHomds.T;
HTHomds.x0 = saddle_point;
HTHomds.extravec = [0 0 0];

HTHomds.UParams = [];    

if dim_npos >= 1
    HTHomds.UParams{1,1} = strcat('UParam',num2str(1));    
    HTHomds.UParams{1,2} = settings.con_UParam1;
end
if dim_npos >= 2
    HTHomds.UParams{2,1} = strcat('UParam',num2str(2));    
    HTHomds.UParams{2,2} = settings.con_UParam2;
end
if dim_npos >= 3
    for f = 3:dim_npos    
        HTHomds.UParams{f,1} = strcat('UParam',num2str(f));        
        HTHomds.UParams{f,2} = 0;        
    end
end

HTHomds.index = 0;


HTHomds.gds.extravec = HTHomds.extravec;
HTHomds.gds.UParams = HTHomds.UParams;
HTHomds.gds.SParams = HTHomds.SParams;




newsettings = settings.copy();
sdata = s(1);
sdata.label = 'HTHom';
sdata.index = 1;
sdata.msg = 'HTHom from simulation';
source = struct('x', x, 'v', [], 'restoreGlobals', @() []);
newsettings.setValue('IP', sdata);
newsettings.getSetting('IP').setSource(source);


newsettings.coord.set(saddle_point);

%  ---> HTHomds -> update into settings .......

%add ntst if not exist
newsettings.addSetting('ntst', CLSetting('ntst', DefaultValues.STARTERDATA.ntst, InputRestrictions.INT_g0, 2, 4, 1, CLSettingsHelp.getHelp('ntst')));
%set ntst
newsettings.setValue('ntst', ntst);

%add ncol if not exists
newsettings.addSetting('ncol', CLSetting('ncol', DefaultValues.STARTERDATA.ncol, InputRestrictions.INT_g0, 2, 4, 2, CLSettingsHelp.getHelp('ncol')));
%set ncol
newsettings.setValue('ncol', ncol);

%add and update ds-structure for homotopy (htds)
newsettings.addSetting('htds', CLSettingHTDS());
newsettings.setValue('htds', HTHomds);
htds = newsettings.getSetting('htds');
htds.cleanupSParam(newsettings);
