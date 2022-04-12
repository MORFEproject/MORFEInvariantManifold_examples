function [newsettings, msg] = GUIselectHTHet(solution)
newsettings = []; msg = '';

settings = solution.settings;
HTHetds = struct();
HTHetds.gds = struct();
saddle_point = solution.connectData.saddle_point;
target_saddle_point = solution.connectData.target_saddle_point;
dim_npos = solution.connectData.dim_npos;
dim = length(saddle_point);
param = solution.param;
[QS, eigvl1] = computeBaseConnecHet(solution.connectData.A,1);%n maal n matrix


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




y1 = zeros(size(solution.y,1),dim+1);
for k = 1:size(solution.y,1)
    y1(k,1:dim) = solution.y(k,:);
    y1(k,dim+1) = norm(target_saddle_point'-solution.y(k,:));
end

x = y1';


D = eig(solution.connectData.A); %dimensie van D : dim maal 1
dim_nneg = sum(real(D) < 0);


ind = size(x,2); %aantal punten
while (ind > 1) && (x(end,ind) - x(end,ind-1) >= 0)
    ind = ind-1;
end

%ind = 0: monotoon stijgend
if ind <= 10
    ind = size(x,2);
end
if ind == size(x,2)
   %monotoon dalend
    verschil = x(end,2:end) - x(end,1:end-1);
    pos = find(verschil > 0);
    if size(pos,2) == 0
        msg = ('No enough data available!');
        return;
    end
end


x=x(1:dim,1:ind);%hier laat je de eps1 weg, hou je enkel de coordinaten over
t=solution.t(1:ind)';
tn=(t-t(1))/(t(end)-t(1));%zo is tn(1) = 0 en tn(end)=1
[a,x,tn] = newmeshHT(x,tn,size(x,2)-1,1,ntst,ncol);
x = interp(tn,1,x,a,ncol);
s=struct();
s(1).data.timemesh = a;
s(1).data.ntst = ntst;
s(1).data.ncol = ncol;
s(1).data.parametervalues = param;
s(1).label = 'ConnecHet';
epsilon1 = norm(x(:,end)-target_saddle_point);

HTHetds.SParams = [];
HTHetds.UParams = [];
for i = 1:(dim-dim_nneg)
    HTHetds.SParams{i,1} = strcat('SParam',num2str(i));    
    HTHetds.SParams{i,2} = 1/epsilon1*(x(:,end)-target_saddle_point)'*QS(:,dim_nneg+i);
end

HTHetds.P0 = param;
HTHetds.ntst = ntst;
HTHetds.ncol = ncol;
HTHetds.nphase = dim;
x = reshape(x,size(x,2)*size(x,1),1);%je maakt er hier 1 lange vector van
HTHetds.T = (t(end)-t(1))/2;%de helft van de periode
HTHetds.eps0 = settings.con_eps0; %hier voeg je epsilon0 toe
HTHetds.eps1 = epsilon1;
s(1).data.T = HTHetds.T;
HTHetds.x0 = saddle_point;
HTHetds.x1 = target_saddle_point;
HTHetds.extravec = [0 0 0];
    
if dim_npos >= 1
    HTHetds.UParams{1,1} = strcat('UParam',num2str(1));    
    HTHetds.UParams{1,2} = settings.con_UParam1;
end
if dim_npos >= 2
    HTHetds.UParams{2,1} = strcat('UParam',num2str(2));    
    HTHetds.UParams{2,2} = settings.con_UParam2;
end
if dim_npos >= 3
    for f = 3:dim_npos    
        HTHetds.UParams{f,1} = strcat('UParam',num2str(f));        
        HTHetds.UParams{f,2} = 0;        
    end
end


HTHetds.TestTolerance = settings.SParamTestTolerance;
         

HTHetds.index = 0;



HTHetds.gds.extravec= HTHetds.extravec;
HTHetds.gds.UParams = HTHetds.UParams;
HTHetds.gds.SParams = HTHetds.SParams;

   
newsettings = settings.copy();
sdata = s(1);
sdata.label = 'HTHet';
sdata.index = 1;
sdata.msg = 'HTHet from simulation';
source = struct('x', x, 'v', [], 'restoreGlobals', @() []);
newsettings.setValue('IP', sdata);
newsettings.getSetting('IP').setSource(source);
newsettings.coord.set(saddle_point);

%  ---> HTHet -> update into settings .......


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
newsettings.setValue('htds', HTHetds);
htds = newsettings.getSetting('htds');
htds.cleanupSParam(newsettings);


