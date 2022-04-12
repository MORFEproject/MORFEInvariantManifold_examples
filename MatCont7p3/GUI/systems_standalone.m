function varargout = systems_standalone(varargin)
% SYSTEM Application M-file for system.fig
%    FIG = SYSTEM launch system GUI.
%    SYSTEM('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 07-Mar-2011 17:47:48
global gds oldgds path_sys MC driver_window;

if nargin == 0 ||((nargin ==1)&&(strcmp(varargin{1},'new'))) % LAUNCH GUI
    h=gcbo;
    arg=get(h,'Tag');
    fig = openfig(mfilename,'reuse', 'invisible');
 % Use system color scheme for figure:
set(fig,'Color',get(0,'DefaultUicontrolBackgroundColor'));

   if strcmp(arg,'systems_standalone')||((nargin ==1)&&strcmp(varargin{1},'new'))
        if ~isempty(MC.starter),gds.open.figuur=1;else gds.open.figuur=0;end
        delete(MC.starter);MC.starter=[];
        % continuer-window open?
        if ~isempty(MC.continuer), gds.open.continuer=1;else gds.open.continuer=0;end
        delete(MC.continuer);MC.continuer=[];
        % numeric window open?   
        if ~isempty(MC.numeric_fig), gds.open.numeric_fig=1;else gds.open.numeric_fig=0;end
        close(MC.numeric_fig);MC.numeric_fig=[];
        %2D-plot open      
        if ~isempty(MC.D2), gds.open.D2=size(MC.D2);else gds.open.D2=0;end
        close(MC.D2);MC.D2=[];
        %3D-plot open      
        if ~isempty(MC.D3), gds.open.D3=size(MC.D3);else gds.open.D3=0;end
        close(MC.D3);MC.D3=[];%   
        %PRC-plot open      
        if ~isempty(MC.PRC), gds.open.PRC=size(MC.PRC);else gds.open.PRC=0;end
        close(MC.PRC);MC.PRC=[];%    
        %dPRC-plot open      
        if ~isempty(MC.dPRC), gds.open.dPRC=size(MC.dPRC);else gds.open.dPRC=0;end
        close(MC.dPRC);MC.dPRC=[];%   
        if ~isempty(MC.integrator),gds.open.integrator=1;else gds.open.integrator=0;end;
        delete(MC.integrator);MC.integrator=[];
        if isfield(gds,'der')
            oldgds=gds;
            init;
        else
            init;
            oldgds=gds;
        end
    end
    % Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
   guidata(fig, handles);
   %checks if the symbolic toolbox is installed
   if (exist('sym')==2)
        set(handles.text14,'String','- symbolically');
        set(handles.f1,'Callback','systems_standalone(''symbolic_callback'',gcbo)');
        set(handles.f2,'Callback','systems_standalone(''symbolic_callback'',gcbo)');
        set(handles.f3,'Callback','systems_standalone(''symbolic_callback'',gcbo)');
        set(handles.f4,'Callback','systems_standalone(''symbolic_callback'',gcbo)');
        set(handles.f5,'Callback','systems_standalone(''symbolic_callback'',gcbo)');       
    end
    load_system(handles);
    if nargout > 0
		varargout{1} = fig;
    end

    % http://undocumentedmatlab.com/blog/customizing-listbox-editbox-scrollbars/
    try % might not work on older versions of matlab
        jScrollPane = findjobj(handles.sys);
        set(jScrollPane, 'HorizontalScrollBarPolicy', ...
            jScrollPane.java.HORIZONTAL_SCROLLBAR_ALWAYS);
        cbStr = sprintf('set(gcbo,''HorizontalScrollBarPolicy'',%d)', ...
            jScrollPane.java.HORIZONTAL_SCROLLBAR_ALWAYS);
        hjScrollPane = handle(jScrollPane,'CallbackProperties');
        set(hjScrollPane,'ComponentResizedCallback',cbStr);
        jViewPort = jScrollPane.getViewport;
        jEditbox = jViewPort.getComponent(0);
        jEditbox.setWrapping(false);
    end

    movegui(fig, 'center');
    fig.Position(3) = fig.Position(3)*1.4;
    set(fig, 'visible', 'on');
    gds.ok = false;
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

    try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch ME
        fprintf(getReport(ME));
        %disp(ME)
		errordlg(ME.message,mfilename);
%         delete(driver_window);  
        global waithndl;
        delete(waithndl);
        waithndl = [];
    end

end

% --------------------------------------------------------------------
function ok_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.ok.
global gds hds path_sys MC;
filterspaces = @(z) z(z ~= ' ');
string_jac='';string_jacp='';string_hess='';string_hessp='';string_tensor3='';string_tensor4='';string_tensor5='';
if isempty(gds.system)
    errordlg('You have to give the system a name','Name');
    return
end



original_cor=filterspaces(get(handles.coordinates,'String'));
original_pa=filterspaces(get(handles.parameters,'String'));
original_equations=deblank(get(handles.sys,'String'));
[equations, cor, pa] = renameforsym(original_equations, original_cor, original_pa);

gds.equations = equations;
gds.parameters = toGdsStruct(pa);
gds.coordinates = toGdsStruct(cor);
par = pa;
if isempty(par)
    button = questdlg('The system has no parameters! Do you want to continue? ',...
    'Parameters','Yes','No','No');
    if strcmp(button,'No')
        return
    end
end
if (~isempty(par))
   par=strcat(',',par);
end
t=get(handles.time,'String');
if (~isempty(t))
    t=strcat(t,',');
end


try
    string_sys=replace_sys_input(gds.equations);
catch
    errordlg('Equations are in the wrong order, compared to the coordinates.','Error')
    return
end
if strcmp(string_sys,'error')
    errordlg('The left-hand sides do not match the coordinates.','Error');
    return
end
fwrite=strcat(gds.system,'.m');
fwrite=fullfile(path_sys,fwrite);
[fid_write,message]=fopen(fwrite,'w');
if fid_write==-1
    errordlg(message,'Error (1)');
    return
end
fread=fullfile(path_sys,'standard.m');

[fid_read,message]=fopen(fread,'r');
if fid_read==-1
    errordlg(message,'Error (2)');
    return
end
global waithndl;
waithndl=waitbar(0,'Precomputing model');
string_handles={'out{1} = @init;';
                'out{2} = @fun_eval;';
                'out{3} = [];';
                'out{4} = [];';
                'out{5} = [];';
                'out{6} = [];';
                'out{7} = [];';
                'out{8} = [];';
                'out{9} = [];';
                'return;';};
string_init=cellstr(make_init);



if (gds.der(3,1)==1||gds.der(4,1)==1)
    str_init=string_init{2};
    str_handles1=string_handles{3};
    str_handles2=string_handles{4};
    str_on='''Jacobian'',handles(3)';str_off='''Jacobian'',[]';
    strp_on='''JacobianP'',handles(4)';strp_off='''JacobianP'',[]';
    if ~isempty(par)
        str_init=strrep(str_init,str_off,str_on);
    end
    str_handles1=strrep(str_handles1,'[]','@jacobian');
    string_handles{3,1}=str_handles1;
    str_handles2=strrep(str_handles2,'[]','@jacobianp');
    string_handles{4,1}=str_handles2;
    str_init=strrep(str_init,strp_off,strp_on);string_init{2,1}=str_init;
    string_jac=gds.jac;string_jacp=gds.jacp;
end
waitbar(0.1);
if (exist('sym')==2&& gds.der(4,1)==1)
    string_jac=cellstr(symjac(handles,gds.equations,cor,pa,10));
    string_jacp=cellstr(symjac(handles,gds.equations,cor,pa,11));
end
waitbar(0.2);
if gds.der(3,1)==1
    string_jac=cellstr(replace_jac_input(gds.jac));
    string_jacp=cellstr(replace_jacp_input(gds.jacp));
end
waitbar(0.3);
if (gds.der(3,2)==1||gds.der(4,2)==1)
    str_init=string_init{2};
    str_handles1=string_handles{5};
    str_handles2=string_handles{6};
    str_on='''Hessians'',handles(5)';str_off='''Hessians'',[]';
    strp_on='''HessiansP'',handles(6)';strp_off='''HessiansP'',[]';
    if ~isempty(par)
    str_init=strrep(str_init,str_off,str_on);
    end
    str_handles1=strrep(str_handles1,'[]','@hessians');
    string_handles{5,1}=str_handles1;
    str_handles2=strrep(str_handles2,'[]','@hessiansp');
    string_handles{6,1}=str_handles2;
    str_init=strrep(str_init,strp_off,strp_on);string_init{2,1}=str_init;
    string_hess=gds.hess;string_hessp=gds.hessp;
end
waitbar(0.4);
if (exist('sym')==2&& gds.der(4,2)==1)
    string_hess=cellstr(symjac(handles,gds.equations,cor,pa,20));
    string_hessp=cellstr(symjac(handles,gds.equations,cor,pa,21));
end
waitbar(0.5);
if (gds.der(3,2)==1)
    string_hess=cellstr(replace_hess_input(gds.hess));
    string_hessp=cellstr(replace_hessp_input(gds.hessp));
end
waitbar(0.6);
if (gds.der(4,3)==1), string_tensor3=gds.tensor3;end
if (exist('sym')==2&& gds.der(4,3)==1)
    string_tensor3=cellstr(symjac(handles,gds.equations,cor,pa,3));
    str_handles2=string_handles{7};
    str_handles2=strrep(str_handles2,'[]','@der3');
    string_handles{7,1}=str_handles2;
end
waitbar(0.7);
if (gds.der(4,4)==1), string_tensor4=gds.tensor4;end
if (exist('sym')==2&& gds.der(4,4)==1)
    string_tensor4=cellstr(symjac(handles,gds.equations,cor,pa,4));
    str_handles2=string_handles{8};
    str_handles2=strrep(str_handles2,'[]','@der4');
    string_handles{8,1}=str_handles2;
end
waitbar(0.8);
if (gds.der(4,5)==1), string_tensor5=gds.tensor5;end
if (exist('sym')==2&& gds.der(4,5)==1)
    string_tensor5=cellstr(symjac(handles,gds.equations,cor,pa,5));
    str_handles2=string_handles{9};
    str_handles2=strrep(str_handles2,'[]','@der5');
    string_handles{9,1}=str_handles2;
end
waitbar(0.9);
if ~isempty(gds.options.UserfunctionsInfo)
    siz = size(gds.options.UserfunctionsInfo,2);
    for i = 1:siz
        string_handles{9+i,1}= sprintf('out{%d}= @%s;',9+i,gds.options.UserfunctionsInfo(i).name);
    end
else siz=0;end

h=0;
filecontent = '';
while feof(fid_read)==0
    tline=fgetl(fid_read);
    h=h+1;
    if h==2
        for i=1:9+siz
            fprintf(fid_write,'%s\n',string_handles{i,1});
            %filecontent = [filecontent,  sprintf('%s\n',string_handles{i,1})];
        end
    end        
    matches=strrep(tline,'time,',t);
    matches=strrep(matches,'odefile',gds.system);
    matches=strrep(matches,',parameters',par);
    fprintf(fid_write,'%s\n',matches); 
    filecontent = [filecontent,  sprintf('%s\n',matches)]; 
    if isfield(gds,'userfunction')
        if ~isempty(findstr(matches,'varargout{1}=der5(coordinates,'))          
            for i = 1:size(gds.userfunction,2)
                hs1 = sprintf('case ''%s''\n\tvarargout{1}=%s(coordinates%s);',gds.options.UserfunctionsInfo(i).name,gds.options.UserfunctionsInfo(i).name,par);
                fprintf(fid_write,'%s\n',hs1); 
                filecontent = [filecontent,  sprintf('%s\n',hs1)]; 
            end
        end
    end        
    if ~isempty(findstr(matches,'function dydt'))
        [dim,x]=size(string_sys);         
        for i=1:dim
              fprintf(fid_write,'%s\n',string_sys{i}); 
              filecontent = [filecontent,  sprintf('%s\n',string_sys{i})];
        end
    end
    if ~isempty(findstr(matches,'handles'))
        [dim,x]=size(string_init);         
        for i=1:dim
              fprintf(fid_write,'%s\n',string_init{i});
              filecontent = [filecontent,  sprintf('%s\n',string_init{i})];
        end
    end
    if (~isempty(findstr(matches,'function jac '))&& ~isempty(string_jac))
        [dim,x]=size(string_jac);         
        for i=1:dim
              fprintf(fid_write,'%s\n',string_jac{i});
              filecontent = [filecontent,  sprintf('%s\n',string_jac{i})];
        end
    end
   
    if (~isempty(findstr(matches,'function jacp'))&& ~isempty(string_jacp))
       [dim,x]=size(string_jacp);       
       for i=1:dim
           fprintf(fid_write,'%s\n',string_jacp{i});
           filecontent = [filecontent,  sprintf('%s\n',string_jacp{i})];
       end
    end
    
    if (~isempty(findstr(matches,'function hess '))&& ~isempty(string_hess))
        [dim,x]=size(string_hess);        
        for i=1:dim
              fprintf(fid_write,'%s\n',string_hess{i});
              filecontent = [filecontent,  sprintf('%s\n',string_hess{i})];
        end
    end
    if (~isempty(findstr(matches,'function hessp'))&& ~isempty(string_hessp))
        [dim,x]=size(string_hessp);
        for i=1:dim
              fprintf(fid_write,'%s\n',string_hessp{i});
              filecontent = [filecontent,  sprintf('%s\n',string_hessp{i})];
        end
    end
    if (~isempty(findstr(matches,'function tens3'))&& ~isempty(string_tensor3))
        [dim,x]=size(string_tensor3);
        for i=1:dim
              fprintf(fid_write,'%s\n',string_tensor3{i});
              filecontent = [filecontent,  sprintf('%s\n',string_tensor3{i})];
        end
    end
    if (~isempty(findstr(matches,'function tens4'))&& ~isempty(string_tensor4))
        [dim,x]=size(string_tensor4);
        for i=1:dim
              fprintf(fid_write,'%s\n',string_tensor4{i});
              filecontent = [filecontent, sprintf('%s\n',string_tensor4{i})];
        end
    end
    if (~isempty(findstr(matches,'function tens5'))&& ~isempty(string_tensor5))
        [dim,x]=size(string_tensor5);
        for i=1:dim
              fprintf(fid_write,'%s\n',string_tensor5{i});
              filecontent = [filecontent, sprintf('%s\n',string_tensor5{i})];
        end
    end
end
newlines = strfind(filecontent, 10);
newline = newlines(1);
gds.filecontent = filecontent(newline+1:end);

if ~isempty(gds.options.UserfunctionsInfo)    
   for i=1:size(gds.options.UserfunctionsInfo,2)
       str_user = []; res=0;
       if isfield(gds,'userfunction') && ~isempty(gds.userfunction{i})
           str_user = systems_standalone('replace_token', cellstr(renameforsym(gds.userfunction{i}, original_cor, original_pa)));
       else 
           str_user=cellstr('res=');
       end
       hs1 = sprintf('function userfun%d=%s(t,kmrgd%s)',i,gds.options.UserfunctionsInfo(i).name,par);
       fprintf(fid_write,'%s\n',hs1);
       hs1 = sprintf('userfun%d',i);
       dim = size(str_user,1);
       for j = 1:dim
           userline = str_user{j};
           d = strmatch('res=',userline,'exact');
           if findstr('res',userline),res=1;end
           userline = strrep(userline,'res',hs1);
           if d==1
               fprintf(fid_write,'\t%s=0;\n',hs1);
           else 
               fprintf(fid_write,'\t%s;\n',userline);
           end
       end
       if res==0,fprintf(fid_write,'\t%s=0;\n',hs1);end
   end
end            
    
waitbar(0.95);
fclose(fid_read);
fclose(fid_write);
file=fullfile(path_sys,gds.system);

%fix mapping here
gds.equations = original_equations;
gds.coordinates = toGdsStruct(original_cor);
gds.parameters = toGdsStruct(original_pa);

save(file,'gds');
gds.ok = true;
hds=[];
delete(waithndl);
delete(handles.system);
%set(MC.mainwindow.compute,'enable','off');
%set(MC.mainwindow.window,'enable','off');
%set(MC.mainwindow.Type,'enable','on');
%set(MC.mainwindow.select_userfunctions,'enable','on');
%{
if gds.open.figuur==1;starter;end
if gds.open.continuer==1;continuer;end
if gds.open.numeric_fig==1;numeric;end
if gds.open.D2>0,for i=1:gds.open.D2,D2;end; end
if gds.open.D3>0,for i=1:gds.open.D3,plotD3;end; end
if gds.open.integrator==1;integrator;end
%}
cd(path_sys);cd ..;
rehash;

% tempstr.label = 'Point';
% tempstr.Tag = 'P_O_DO';
% matcont('point_callback',tempstr)

% --------------------------------------------------------------------
function varargout = cancel_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.cancel.
global gds oldgds path_sys;
gds=oldgds;
if (~isempty(gds))
    file=fullfile(path_sys,gds.system);
    save(file,'gds');
    load_system(handles);
    %if gds.open.figuur==1;starter;end
    %if gds.open.continuer==1;continuer;end
    %if gds.open.numeric_fig==1;numeric;end
    %if gds.open.D2>0,for i=1:gds.open.D2,D2;end; end
    %if gds.open.D3>0,for i=1:gds.open.D3,plotD3;end; end
    %if gds.open.integrator==1;integrator;end
end
delete(handles.system);

% --------------------------------------------------------------------
function varargout = coordinates_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.coordinates.
global gds;
% store old coordinates, they might contain usefull initial values
old_coord = gds.coordinates;
gds.coordinates=[];gds.plot2='';gds.plot3='';gds.PRC='';gds.dPRC='';
str=get(handles.coordinates,'String');
str=parse(str);
t=1;
for i=1:length(str)
    co{t,1}=str{i,1};
    string=str{i};
    str1=findstr(string,'[');
    str2=findstr(string,']');
    m=length(str1);
    if (m==length(str2)&&(m==1))
        num=str2double(string(str1+1:str2-1));
        var=string(1:str1-1);
        co{t,1}=strcat(var,'(1)');
        if num>1
            for j=2:num
                t=t+1;
                co{t,1}=strcat(var,'(',num2str(j),')');
            end
        end
    end
    t=t+1;
end
gds.coordinates=co;
gds.dim=size(gds.coordinates,1);
if ((gds.dim==1)&&(strcmp(gds.coordinates{1},'')))
    gds.coordinates=[];
    gds.dim=0;
else
    % run through each coordinate and try to set its initial value
    % from the old coordinate
    for i=1:gds.dim
        gds.coordinates{i,2}=0;
        for j = 1:length(old_coord)
            if strcmp(old_coord{j,1},gds.coordinates{i,1})
                gds.coordinates{i,2}=old_coord{j,2};
                break;
            end
        end
    end
end




% --------------------------------------------------------------------
function varargout = parameters_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.parameters.
global gds;
% store old parameters, they might contain usefull initial values
old_params = gds.parameters;
gds.parameters=[];gds.options.ActiveParams=[];
gds.T=[];gds.eps0=[];gds.eps1=[];gds.extravec=[1 0 0]; gds.t=[]; gds.epsilon=[];
gds.options.IgnoreSingularity=[];gds.plot2='';gds.plot3='';gds.PRC='';gds.dPRC='';
str=get(handles.parameters,'String');
gds.parameters=parse(str);
[ndim,x]=size(gds.parameters);


% display a warning to the user to also change relevant jacobians..
% but only if th user has set them by hand, to check this, we look into
% the gds.der matrix, if the relevant fields have been set, the find
% will return them and so not be empty
if (~isempty(find(gds.der(2:3,1:5),1))&&(~isempty(gds.jacp)||(~isempty(gds.hessp))))
    str=sprintf('Don''t forget to change jacobianp/hessianp \n you also have to change ''df..d..='' etc.');
    warndlg(str,'warning!!');
end
if ((ndim==1)&&(strcmp(gds.parameters{1},'')))
    gds.parameters=[];
else
    % run through each parameter and try to set its initial value
    % from the old parameter
    for i=1:ndim
        gds.parameters{i,2}=0;
        for j = 1:length(old_params)
            if strcmp(old_params{j,1},gds.parameters{i,1})
                gds.parameters{i,2}=old_params{j,2};
                break;
            end
        end
    end
end

  
% --------------------------------------------------------------------
function varargout = time_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.time.
global gds
gds.time=[];
str=get(handles.time,'String');
gds.time={str, 0};

% --------------------------------------------------------------------
function varargout = name_system_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.name_system.
global gds path_sys
str1=get(handles.name_system,'String');
str1=strrep(str1,'.m','');
str=strcat(str1,'.m');
file=fullfile(path_sys,str);
if (exist(file)==2)
    warning_box(str1,file,handles);
else 
    w=which(str);
    if ~isempty(w)
        errordlg('This name is already in use for another function (Possibly in Matlab itself).','wrong name');
        set(handles.name_system,'String','');gds.system='';
        return
    end
    gds.system=str1;
end

% --------------------------------------------------------------------
function varargout = edit8_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit8
global gds
gds.equations=[];

% --------------------------------------------------------------------
function varargout = numerically_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.n1.
global gds
t=get(h,'Tag');
num=str2double(strtok(t,'n'));
set(h,'Value',1);
gds.der(1,num)=1;gds.der(2,num)=0;
gds.der(3,num)=0;gds.der(4,num)=0;
set(0,'ShowHiddenHandles','on');
hf=findobj('Tag',strcat('f',num2str(num)));
hr=findobj('Tag',strcat('r',num2str(num)));
val=get(hr,'Value');
if (val==1)
  edit=findobj('Tag','edit');  stat=findobj('Tag','stat');
  editp=findobj('Tag','editp');  statp=findobj('Tag','statp');
  delete(edit);  delete(stat);
  delete(editp);  delete(statp);
end
set(0,'ShowHiddenHandles','off');
pos=[0.0335 0.09 0.925 0.47];
set(handles.sys,'Position',pos);
set(hf,'Value',0);set(hr,'Value',0);
% --------------------------------------------------------------------
function varargout = routine_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.r1.
global gds
p=0;jac='';jacp='';hess='';hessp='';
for j=1:gds.dim
    for i=1:gds.dim
        p=p+1;
        jac{p,1}=strcat('df',gds.coordinates{j,1},'d',gds.coordinates{i,1},'=');
    end
end
p=0;
for j=1:gds.dim
    for i=1:gds.dim
        for m=1:gds.dim
            p=p+1;
            hess{p,1}=strcat('df',gds.coordinates{j,1},'d',gds.coordinates{i,1},'d',gds.coordinates{m,1},'=');
        end
    end     
end
dimp=size(gds.parameters,1);
if (dimp==0)
    jacp='';hessp='';
else
    p=0;
    for j=1:gds.dim
        for i=1:dimp
            p=p+1;
            jacp{p,1}=strcat('df',gds.coordinates{j,1},'d',gds.parameters{i,1},'=');
        end
    end
    p=0;
    for j=1:gds.dim
        for i=1:gds.dim
            for m=1:dimp
                p=p+1;
                hessp{p,1}=strcat('df',gds.coordinates{j,1},'d',gds.coordinates{i,1},'d',gds.parameters{m,1},'=');
            end
        end     
    end
end
set(h,'Value',1);
t=get(h,'Tag');
num=str2double(strtok(t,'r'));
set(0,'ShowHiddenHandles','on');
hf=findobj('Tag',strcat('f',num2str(num)));
hn=findobj('Tag',strcat('n',num2str(num)));
hn=findobj('Tag',strcat('n',num2str(num)));
edit=findobj('Tag','edit');stat=findobj('Tag','stat');
editp=findobj('Tag','editp');statp=findobj('Tag','statp');
delete(edit);delete(stat);
delete(editp);delete(statp);
set(0,'ShowHiddenHandles','off');
gds.der(3,num)=1;gds.der(1,num)=0;
gds.der(2,num)=0;gds.der(4,num)=0;
set(hf,'Value',0);set(hn,'Value',0);
pos=[0.0335 0.34 0.925 0.23];
set(handles.sys,'Position',pos);
color=get(0,'defaultUicontrolBackgroundColor');
switch num
        case 1, 
            string='1st order derivatives'; stringp='Enter the jacobianp';
        case 2,
             string='2nd order derivatives';stringp='Enter the hessianp';
        case 3,
             string='3thd order derivatives';
        otherwise, 
             string='4th order derivatives';    
end
poss=  [0.0335 0.31 0.925 0.029];pose=  [0.0335 0.09 0.925 0.22];
posje1=[0.0335 0.21 0.925 0.10];posjs= [0.0335 0.18 0.925 0.03];
posje2=[0.0335 0.09 0.925 0.09];
stat=uicontrol(handles.system,'Style','text','HorizontalAlignment','left','String',string,'Tag','stat','BackGroundColor',color,'Units','normalized','Position',poss,'fontname','FixedWidth','fontsize',12);
%if ((num==1)&&~isempty(gds.jac))
 %       str=gds.jac; strp=gds.jacp;
 %els
 if (num==1)
    gds.jac='';gds.jacp='';
    str=jac; strp=jacp;
        %elseif ((num==2)&& ~isempty(gds.hess))
       % str=gds.hess; strp=gds.hessp;
elseif (num==2)
   gds.hess='';gds.hessp='';
   str=hess; strp=hessp;
else
   str=''; strp='';
end
edit=uicontrol(handles.system,'Style','edit','String',str,'HorizontalAlignment','left','Tag','edit','BackGroundColor',[1 1 1],'Max',40,'Callback','der_callback','Userdata',num,'Units','normalized','Position',pose,'fontname','FixedWidth','fontsize',12);
if ((num==1)&&(dimp~=0))||((num==2)&&(dimp~=0))
  edit=uicontrol(handles.system,'Style','edit','String',str,'HorizontalAlignment','left','Tag','edit','BackGroundColor',[1 1 1],'Max',40,'Callback','der_callback','Userdata',num,'Units','normalized','Position',posje1,'fontname','FixedWidth','fontsize',12); 
  stat_p=uicontrol(handles.system,'Style','text','HorizontalAlignment','left','String',stringp,'Tag','statp','BackGroundColor',color,'Units','normalized','Position',posjs);  
  edit_p=uicontrol(handles.system,'Style','edit','String',strp,'HorizontalAlignment','left','Tag','editp','BackGroundColor',[1 1 1],'Max',40,'Callback','der_callback','Userdata',num,'Units','normalized','Position',posje2,'fontname','FixedWidth','fontsize',12);
end

%-------------------------------------------------------------------
function varargout = file_callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.r1.
global gds
p=0;
t=get(h,'Tag');
num=str2double(strtok(t,'f'));
switch num
case 1
    prompt  = {'Enter the file name of the jacobian:','Enter the file name of the jacobianp:'};
    title   = 'Input name of derivative-file';
    lines= 1;
    def     = {'',''};
    answer  = inputdlg(prompt,title,lines,def);
    if isempty(answer)||isempty(answer{1})||isempty(answer{2})
       set(h,'Value',0);
       return;
    end
    jacobian=answer{1};jacobianp=answer{2};
    gds.jac=parsefile(jacobian);
    gds.jacp=parsefile(jacobianp);
    gds.jac=strcat('jac=',replace_token(gds.jac),';');gds.jacp=strcat('jacp=',replace_token(gds.jacp),';');
case 2
    prompt  = {'Enter the file name of the hessian:','Enter the file name of the hessianp:'};
    title   = 'Input name of derivative-file';
    lines= 1;
    def     = {'',''};
    answer  = inputdlg(prompt,title,lines,def);
    if isempty(answer)|| isempty(answer{1})||isempty(answer{2})
       set(h,'Value',0);
       return;
    end
    hessian=answer{1};hessianp=answer{2};
    gds.hess=[];gds.hessp=[];
    hess=parsefile(hessian);hess=replace_token(hess);
    hessp=parsefile(hessianp);hessp=replace_token(hessp);
    for i=1:gds.dim
        gds.hess{i,1}=strcat('hess',num2str(i),'=',hess{i,1},';');
        gds.hess{gds.dim+i,1}=strcat('hess','(:,:,',num2str(i),') =','hess',num2str(i),';');
    end
    dimp=size(gds.parameters,1);
    for  i=1:dimp
        gds.hessp{i,1}=strcat('hessp',num2str(i),'=',hessp{i,1},';');
        gds.hessp{dimp+i,1}=strcat('hessp','(:,:,',num2str(i),') =','hessp',num2str(i),';');
    end
case {3,4,5}
    prompt  = {'Enter the file name of the tensor:'};
    title   = 'Input name of derivative-file';
    lines= 1;
    def     = {''};
    answer  = inputdlg(prompt,title,lines,def);
    if isempty(answer)||isempty(answer{1})
       set(h,'Value',0);
       return;
    end
    tensor=answer{1};
    tens=parsefile(tensor);tens=replace_token(tens);
    switch num
    case 3
        dim=gds.dim*gds.dim;
        gds.tensor3=[];j=0;
        for i=1:dim
            if (mod(i,gds.dim)==1),j=j+1;end
            k=mod(i,gds.dim);if (k==0),k=gds.dim;end
            gds.tensor3{i,1}=strcat('tens3',num2str(i),'=',tens{i},';');
            gds.tensor3{dim+i,1}=strcat('tens3','(:,:,',num2str(j),',',num2str(k),') =','tens3',num2str(i),';');
        end
    case 4
        dim=(gds.dim)^3;
        gds.tensor4=[];j=0;d=0;
        for  i=1:dim
            k=mod(i,gds.dim);if (k==0),k=gds.dim;end
            if (mod(i,gds.dim)==1),j=j+1;end
            j=mod(j,gds.dim);if j==0,j=gds.dim;end
            if (mod(i,gds.dim*gds.dim)==1),d=d+1;end
            gds.tensor4{i,1}=strcat('tens4',num2str(i),'=',tens{i},';');
            gds.tensor4{dim+i,1}=strcat('tens4','(:,:,',num2str(d),',',num2str(j),',',num2str(k),') =','tens4',num2str(i),';');
        end
    case 5
        dim=(gds.dim)^4;
        gds.tensor4=[];j=0;d=0;p=0;
        for i=1:dim
            k=mod(i,gds.dim);if (k==0),k=gds.dim;end
            if (mod(i,gds.dim)==1),j=j+1;end
            j=mod(j,gds.dim);if j==0,j=gds.dim;end
            if (mod(i,gds.dim*gds.dim)==1),d=d+1;end
            d=mod(d,gds.dim);if d==0,d=gds.dim;end
            if (mod(i,(gds.dim)^3)==1),p=p+1;end
            gds.tensor5{i,1}=strcat('tens5',num2str(i),'=',tens{i},';');
            gds.tensor4{dim+i,1}=strcat('tens5','(:,:,',num2str(p),',',num2str(d),',',num2str(j),',',num2str(k),') =','tens5',num2str(i),';');
        end
    end
end
jac='';jacp='';hess='';hessp='';tens='';
set(h,'Value',1);t=get(h,'Tag');
set(0,'ShowHiddenHandles','on');
hr=findobj('Tag',strcat('r',num2str(num)));
hn=findobj('Tag',strcat('n',num2str(num)));
val=get(hr,'Value');
if (val==1)   
    edit=findobj('Tag','edit');stat=findobj('Tag','stat');
    editp=findobj('Tag','editp');statp=findobj('Tag','statp');
    delete(edit);delete(stat);delete(editp);delete(statp);
end
set(0,'ShowHiddenHandles','off');
pos=[0.0335 0.09 0.925 0.47];
set(handles.sys,'Position',pos);
gds.der(4,num)=1;gds.der(1,num)=0;
gds.der(2,num)=0;gds.der(3,num)=0;
set(hr,'Value',0);set(hn,'Value',0);  
    
%-------------------------------------------------------------------
function varargout = symbolic_callback(h)
% Stub for Callback of the uicontrol handles.r1.
global gds
p=0;jac='';jacp='';hess='';hessp='';tens='';
set(h,'Value',1);t=get(h,'Tag');
num=str2double(strtok(t,'f'));
set(0,'ShowHiddenHandles','on');
hr=findobj('Tag',strcat('r',num2str(num)));
hn=findobj('Tag',strcat('n',num2str(num)));
hsys=findobj('Tag','sys');
val=get(hr,'Value');
if (~isempty(val) && (val==1))   
    edit=findobj('Tag','edit');stat=findobj('Tag','stat');
    editp=findobj('Tag','editp');statp=findobj('Tag','statp');
    delete(edit);delete(stat);delete(editp);delete(statp);
end
set(0,'ShowHiddenHandles','off');
pos=[0.0335 0.09 0.925 0.47];
set(hsys,'Position',pos);
gds.der(4,num)=1;gds.der(1,num)=0;
gds.der(2,num)=0;gds.der(3,num)=0;
set(hr,'Value',0);set(hn,'Value',0);

%-------------------------------------------------------------------
function sym_jac = symjac(handles,string,cor,p1e,number)
% Stub for Callback of the uicontrol handles.r1.
global gds;
%string=get(handles.sys,'String');
dimp=size(gds.parameters,1);
if isempty(string)
    string = '';
    return
end
string_sys = cellstr(string);
string=''; tempory = '';equation = '';
for j1e = 1:(gds.dim)
    syss{j1e} = strcat(gds.coordinates{j1e,1},char(39));
end
[tempory,equation] = parse_input(string_sys,syss);
num_temp = size(tempory,1);
for i1e = 1:num_temp
    string{i1e,1} = strcat(tempory{i1e,1},';');
end
for i = 1:num_temp
    assignment = tempory{i,1};
    equation_sign = strfind(assignment,'=');
    lhs_assignment = strtrim(assignment(1:equation_sign-1));
    syms(lhs_assignment);
end
num_eq = size(equation,1);
for i1e=1:num_eq
    string{num_temp+i1e,1} = '';
end
for j1e = 1:num_eq
      [null,string{num_temp+j1e,1}] = strtok(equation{j1e,1},'=');
      string{num_temp+j1e,1} = strtok(string{num_temp+j1e,1},'=');
      string{num_temp+j1e,1} = strcat(string{num_temp+j1e,1},';');
end
string{num_temp+1,1} = strcat('dydt=[',string{num_temp+1,1});
string{num_temp+num_eq,1} = strcat(string{num_temp+num_eq,1},'];');
if (num_eq~=gds.dim)
  string='error';
  return;
end
strings=char(string);
cor=strrep(cor,',',' ');
cors=sprintf('syms %s',cor);eval(cors);
p1a=strrep(p1e,',',' ');
p1e=sprintf('syms %s',p1a);eval(p1e);
stri='';
for i1ee=1:size(strings,1)
    stri=strcat(stri,strings(i1ee,:));
end
eval(stri);
cors=strcat('[',cor,']');
jac=jacobian(dydt,eval(cors));
if number==10% case jacobian

    j1c = symmat2line(jac);
    sym_jac=replace_token(cellstr(strcat('jac=',j1c,';')));
elseif number==11 %case jacobianp
    p1a=strcat('[',p1a,']');
    j1p= symmat2line( jacobian(dydt,eval(p1a)) );    
    sym_jac=replace_token(cellstr(strcat('jacp=',j1p,';')));
   
elseif number==20 %case hessians
    gds.hess=[];
    for i1e=1:gds.dim
        h1s=symmat2line( diff(jac,eval(gds.coordinates{i1e,1})) );
        gds.hess{i1e,1}=strcat('hess',num2str(i1e),'=',h1s,';');
        gds.hess{gds.dim+i1e,1}=strcat('hess','(:,:,',num2str(i1e),') =','hess',num2str(i1e),';');
    end
    sym_jac=replace_token(cellstr(gds.hess));
elseif number==21 % case hessianp
    gds.hessp=[];
    for  i1e=1:dimp
        h1s=diff(jac,eval(gds.parameters{i1e,1}));
        h1sp = symmat2line(   diff(jac,eval(gds.parameters{i1e,1}))    );
        gds.hessp{i1e,1}=strcat('hessp',num2str(i1e),'=',h1sp,';');
        gds.hessp{dimp+i1e,1}=strcat('hessp','(:,:,',num2str(i1e),') =','hessp',num2str(i1e),';');
    end
    sym_jac=replace_token(cellstr(gds.hessp));
    
elseif number==3
    dim=gds.dim*gds.dim;
    gds.tensor3=[];je1=0;if(dim==1),je1=1;end
    for  i1e=1:dim
        if (mod(i1e,gds.dim)==1),je1=je1+1;end
        ke1=mod(i1e,gds.dim);
        if (ke1==0),ke1=gds.dim;end
        he1=diff(jac,eval(gds.coordinates{je1,1}));
        h2s=diff(he1,eval(gds.coordinates{ke1,1}));
	h2s = symmat2line( h2s );
        gds.tensor3{i1e,1}=strcat('tens3',num2str(i1e),'=',h2s,';');
        gds.tensor3{dim+i1e,1}=strcat('tens3','(:,:,',num2str(je1),',',num2str(ke1),') =','tens3',num2str(i1e),';');
    end 
    sym_jac=replace_token(cellstr(gds.tensor3));
elseif number==4
    dim=(gds.dim)^3;
    gds.tensor4=[];je1=0;de1=0;if(dim==1),je1=1;de1=1;end
    for  i1e=1:dim
        ke1=mod(i1e,gds.dim);if (ke1==0),ke1=gds.dim;end
        if (mod(i1e,gds.dim)==1),je1=je1+1;end
        je1=mod(je1,gds.dim);if je1==0,je1=gds.dim;end
        if (mod(i1e,gds.dim*gds.dim)==1),de1=de1+1;end
        he1=diff(jac,eval(gds.coordinates{de1,1}));
        h2s=diff(he1,eval(gds.coordinates{je1,1}));
        h3s=diff(h2s,eval(gds.coordinates{ke1,1}));

	h3s = symmat2line( h3s );

        gds.tensor4{i1e,1}=strcat('tens4',num2str(i1e),'=',h3s,';');
        gds.tensor4{dim+i1e,1}=strcat('tens4','(:,:,',num2str(de1),',',num2str(je1),',',num2str(ke1),') =','tens4',num2str(i1e),';');
    end 
    sym_jac=replace_token(cellstr(gds.tensor4));
elseif number==5
    dim=(gds.dim)^4;
    gds.tensor5=[];je1=0;de1=0;pe1=0;if(dim==1),je1=1;de1=1;pe1=1;end
    for  i1e=1:dim
        ke1=mod(i1e,gds.dim);if (ke1==0),ke1=gds.dim;end
        if (mod(i1e,gds.dim)==1),je1=je1+1;end
        je1=mod(je1,gds.dim);if je1==0,je1=gds.dim;end
        if (mod(i1e,gds.dim*gds.dim)==1),de1=de1+1;end
        de1=mod(de1,gds.dim);if de1==0,de1=gds.dim;end
        if (mod(i1e,gds.dim*gds.dim*gds.dim)==1),pe1=pe1+1;end
        he1=diff(jac,eval(gds.coordinates{pe1,1}));
        h2s=diff(he1,eval(gds.coordinates{de1,1}));
        h3s=diff(h2s,eval(gds.coordinates{je1,1}));
        h4s=diff(h3s,eval(gds.coordinates{ke1,1}));
	
	h4s = symmat2line( h4s );

        gds.tensor5{i1e,1}=strcat('tens5',num2str(i1e),'=',h4s,';');
        gds.tensor5{dim+i1e,1}=strcat('tens5','(:,:,',num2str(pe1),',',num2str(de1),',',num2str(je1),',',num2str(ke1),') =','tens5',num2str(i1e),';');
    end 
    sym_jac=replace_token(cellstr(gds.tensor5));
end
    

%--------------------------------------------------------------------
function varargout=load_system(handles)
global gds
co='';par='';
if ((gds.dim)~=0)
    co=gds.coordinates{1,1};
end
if ((gds.dim)>=2)
  for i=2:gds.dim
       co=horzcat(co,',',gds.coordinates{i,1});
  end
end
dim=size(gds.parameters,1);
if (dim~=0)
    par=gds.parameters{1,1};
end
if (dim>=2)
   for i=2:dim
        par=horzcat(par,',',gds.parameters{i,1});
   end
end  
set(handles.name_system,'String',gds.system);
set(handles.coordinates,'String',co);
set(handles.parameters,'String',par);
set(handles.time,'String',gds.time{1,1});
% before setting the system, trim the equations
equations=gds.equations;
equations_string = '';
for equation = equations'
    equations_string = sprintf('%s%s\n',equations_string, strtrim(equation'));
end
set(handles.sys,'String',equations_string);
set(handles.n1,'Value',gds.der(1,1));
set(handles.n2,'Value',gds.der(1,2));
set(handles.n3,'Value',gds.der(1,3));
set(handles.n4,'Value',gds.der(1,4));
set(handles.n5,'Value',gds.der(1,5));
set(handles.f1,'Value',gds.der(4,1));
set(handles.f2,'Value',gds.der(4,2));
set(handles.f3,'Value',gds.der(4,3));
set(handles.f4,'Value',gds.der(4,4));
set(handles.f5,'Value',gds.der(4,5));
set(handles.r1,'Value',gds.der(3,1));
if (gds.der(3,1)==1)
    routine_Callback(handles.r1,[],guidata(handles.r1));
end
set(handles.r2,'Value',gds.der(3,2));
if (gds.der(3,2)==1)
    routine_Callback(handles.r2,[],guidata(handles.r2));
end
guidata(handles.system,handles);

%---------------------------------------------------------------------
function all=parse(input_string)
global gds;
remainder=input_string;
all='';
while (any(remainder))
[chopped,remainder]=strtok(remainder,',');
all=strvcat(all,chopped);
end
all=cellstr(all);

%--------------------------------------------------------------------
function string=replace_sys_input(string)
global gds;
if isempty(string)
    string = '';
    return
end
string_sys = cellstr(string);
string=''; temp = '';eq = '';
for j = 1:(gds.dim)
    sys{j} = strcat(gds.coordinates{j,1},char(39));
end
[temp,eq] = parse_input(string_sys,sys);
num_temp = size(temp,1);
for i = 1:num_temp
    string{i,1} = strcat(temp{i,1},';');
end
num_eq = size(eq,1);
for i=1:num_eq
    string{num_temp+i,1} = '';
end
for j = 1:num_eq
      [null,string{num_temp+j,1}] = strtok(eq{j,1},'=');
      string{num_temp+j,1} = strtok(string{num_temp+j,1},'=');
      string{num_temp+j,1} = strcat(string{num_temp+j,1},';');
end
string{num_temp+1,1} = strcat('dydt=[',string{num_temp+1,1});
string{num_temp+num_eq,1} = strcat(string{num_temp+num_eq,1},'];');
if (num_eq~=gds.dim)
  string='error';
  return;
end
string = replace_token(string);

%--------------------------------------------------------------------
function string = replace_jac_input(string)
global gds;
if isempty(string)
    string = '';
    return
end
string_jac = cellstr(string);
string=''; temp = ''; eq = '';
p = 0;
for j = 1:gds.dim
    for i = 1:gds.dim
        p = p+1;
        jac{p,1} = strcat('df',gds.coordinates{j,1},'d',gds.coordinates{i,1});
    end
end
[temp,eq] = parse_input(string_jac,jac);
string = make_mat(temp,eq,jac,gds.dim,'jac');
string = replace_token(string);

%--------------------------------------------------------------------
function string = replace_jacp_input(string)
global gds;
if isempty(string)
    string = '';
    return
end
string_jacp = cellstr(string);
string=''; temp = ''; eq ='';
p = 0;
dimp = size(gds.parameters,1);
for j = 1:gds.dim
    for i = 1:dimp
        p = p+1;
        jacp{p,1} = strcat('df',gds.coordinates{j,1},'d',gds.parameters{i,1});
    end
end
[temp,eq] = parse_input(string_jacp,jacp);
string = make_mat(temp,eq,jacp,dimp,'jacp');
string = replace_token(string);

%-----------------------------------------------------------------------------
function string = replace_hess_input(string)
global gds;
if isempty(string)
    string = '';
    return
end
string_hess = cellstr(string);
string ='';temp =''; eq ='';
p = 0;
for j = 1:gds.dim
    for i = 1:gds.dim
        for m = 1:gds.dim
            p = p+1;
            hess{p,1} = strcat('df',gds.coordinates{j,1},'d',gds.coordinates{i,1},'d',gds.coordinates{m,1});
        end
    end     
end
[temp,eq] = parse_input(string_hess,hess);
string = make_mathess(temp,eq,hess,gds.dim,'hess');
string = replace_token(string);

%-----------------------------------------------------------------------
function string=replace_hessp_input(string)
global gds;
if isempty(string)
    string = '';
    return
end
string_hessp = cellstr(string);
string=''; temp =''; eq ='';
dimp = size(gds.parameters,1);
p = 0;
for j = 1:gds.dim
    for i = 1:gds.dim
        for m = 1:dimp
            p = p+1;
            hessp{p,1} = strcat('df',gds.coordinates{j,1},'d',gds.coordinates{i,1},'d',gds.parameters{m,1});
        end
    end     
end
[temp,eq] = parse_input(string_hessp,hessp);
string = make_mathess(temp,eq,hessp,dimp,'hessp');
string = replace_token(string);

%--------------------------------------------------------------------
function string = replace_token(string)
global gds
token = strcat('[],;=/*-+ ^()!"',char(39));
dim = size(string,1);
for t=1:dim
    for i=1:(gds.dim)
        repl = string{t,1};
        m = findstr(repl,gds.coordinates{i,1});
        p = length(m);
        h = length(gds.coordinates{i,1});
        j = num2str(i);
        x = strcat('kmrgd(',j,')');
        len = length(x);
        w = 0;
        if p>0
            for g = 1:p
                if (m(g) > 1)  
                    k = findstr(repl(m(g)-1+w*(len-h)),token);
                else
                    k = 1;
                end
                if ((m(g)+h+w*(len-h)) < length(repl))                 
                    s = findstr(repl(m(g)+h+w*(len-h)),token);
                else
                    s = 1;
                end
                if (~isempty(s) && ~isempty(k))
                    repl((m(g)+w*(len-h)):(m(g)+h-1+w*(len-h))) = '';
                    if(m(g)>1)
                        y = m(g)-1+w*(len-h);
                        repl = sprintf('%s%s%s',repl(1:(m(g)-1+w*(len-h))),x,repl(m(g)+w*(len-h):end));
                    else 
                        repl = sprintf('%s%s',x,repl(m(g)+w*(len-h):end));
                    end
                    string{t,1} = repl;
                    w = w+1;
                end
            end
        end
    end
end
%---------------------------------------------------------------------------
function [temp,eq] = parse_input(string,type)
% can't parse empty strings, so first remove them here
j=1;
cleaned_string={};
for i=1:size(string,1)
    if ~strcmp(string{i},'')
        cleaned_string{j}=string{i};
        j=j+1;
    end
end
string=cellstr(cleaned_string');
% continue with code
dim = size(string,1);
vars = size(type,2);
temp='';eq='';
p=1;s=1;
for j=1:dim
    k=[];  
    for i=1:length(type)
        teststring = string{j};
        if exist('strtrim','builtin')
            coordinate = strtrim(type{i});
        else
            coordinate = deblank(type{i});
        end
        match = strcat('\<',coordinate,'\>');
        pos = regexp(teststring,match);
        if ~isempty(pos) && pos(1)==1
           k = 1;
           tmpv = type{1};
           if (vars-i ~= dim-j) && (tmpv(end) == '''')
               error('Equations are in the wrong order, compared to the coordinates.');
           end
        end
    end
    if (findstr(string{j},'=')&(isempty(k)))
        temp{p,1} = string{j};
        h = 1; c = 0; p = p+1;
     elseif (findstr(string{j},'=')&~(isempty(k)))
             eq{s,1} = string{j};
             h = 0; c = 1; s = s+1;
     elseif (~findstr(string{j},'='))
          if (h==1)
             temp{p,1} = strcat(temp{p,1},string{j});
             p = p+1;
          elseif (c==1)
                  eq{s,1} = strcat(eq{s,1},string{j});
                  s = s+1;
          end
    end
end

%-------------------------------------------------------------------------
function string=make_mat(temp,eq,type,dim,str)
global gds;
num_temp=size(temp,1);
for i=1:num_temp
    string{i,1} = strcat(temp{i,1},';');
end
for i=1:(gds.dim)
    string{num_temp+i,1} = '';
end
num_eq = size(eq,1);
t = 0 ; p = 1;
for j=1:length(type)
    for i=1:num_eq        
         x = findstr(eq{i,1},type{j,1});
         if ~isempty(x)
            [null,eq{i,1}] = strtok(eq{i,1},'=');
            eq{i,1} = strtok(eq{i,1},'=');
            if isempty(eq{i,1})
               eq{i,1} = 0;
            end
            t = t+1;
            if (t==(dim+1))
                p = p+1;
                t = 1;
            end
            string{num_temp+p,1} = sprintf('%s%c%c%c%s',string{num_temp+p,1},char(32),char(32),char(32),eq{i,1});
            eq{i,1} = '';    
        end 
    end
end        
for i=1:gds.dim
    string{num_temp+i,1} = strcat(string{num_temp+i},';');
end
string{num_temp+1,1} = strcat(str,'=[',string{num_temp+1,1});
string{num_temp+(gds.dim),1} = strcat(string{num_temp+(gds.dim),1},'];');

%-----------------------------------------------------------------------
function string=make_mathess(temp,eq,type,dim,str)
global gds
num_temp = size(temp,1);
if num_temp>0
   for i=1:num_temp
       string{i} = strcat(temp{i,1},';');
   end
end
for i=1:(dim*gds.dim+dim)
    string{num_temp+i,1} ='';
end
num_eq = size(eq,1);
eq = sort(eq);
p = 1; r = 0; t = 0;
for m=1:dim
    for j=m:dim:length(type)
        for i=1:num_eq
            x = findstr(eq{i,1},type{j,1});
            if ~isempty(x)
               [null,eq{i,1}] = strtok(eq{i,1},'=');    
               eq{i,1} = strtok(eq{i,1},'=');
               if isempty(eq{i,1})
                  eq{i,1} = 0;
               end
               t = t+1;
               if (t==(gds.dim+1))
                  p = p+1;
                  t = 1;
               end
               string{num_temp+p,1} = sprintf('%s%c%c%s',string{num_temp+p,1},char(32),char(32),eq{i,1});
               eq{i,1} = '';
           end 
        end
    end        
    for i=1:gds.dim
        r=r+1;  
        string{num_temp+r,1} = strcat(string{num_temp+r},';');
    end
    string{num_temp+r-(gds.dim)+1,1} = strcat(str,num2str(m),'=[',string{num_temp+r-(gds.dim)+1,1});
    string{num_temp+r,1} = strcat(string{num_temp+r,1},'];');
end
for i=1:dim
string{num_temp+r+i,1} = strcat(str,'(:,:,',num2str(i),') =',str,num2str(i),';');
end

%------------------------------------------------------------------------
function init
global gds;
    gds = []; gds.coordinates = []; gds.parameters = [];
    gds.time{1,1} = 't';gds.time{1,2} = 0; gds.options = contset;
    gds.system = '';
    gds.curve.new = '';gds.curve.old = '';
    gds.equations = [];
    gds.dim = 0;
    gds.der = [[1 1 1 1 1];zeros(3,5)]; 
    gds.jac = '';%string that contains the jacobian
    gds.jacp = '';%string that contains the jacobianp
    gds.hess = '';%string that contains the hessian
    gds.hessp = '';%string that contains the hessianp
    gds.tensor3 = ''; gds.tensor4 = ''; gds.tensor5 = '';
    gds.point = ''; gds.type = '';
    gds.discretization.ntst = 20; gds.discretization.ncol = 4;
    gds.period = 1;
    gds.plot2 = '';gds.plot3 = '';gds.PRC='';gds.dPRC='';
    gds.open.figuur = 0; gds.open.continuer = 0; gds.open.numeric_fig = 0;
    gds.open.D2 = 0;gds.open.D3 = 0;gds.open.PRC = 0; gds.open.dPRC = 0; gds.open.integrator = 0;
    gds.integrator = []; gds.integrator.method = 'ode45'; gds.integrator.options = [];
    gds.integrator.tspan = [0 1]; gds.numeric = [];
    gds.numeric.O = {'time' 1;'coordinates' 1;'parameters' 0'};
    gds.numeric.EP = {'coordinates' 1;'parameters' 1;'testfunctions' 0;'eigenvalues' 0;'current stepsize' 0};
    % XXXX
    %gds.numeric.LC = {'parameters' 1;'period' 1;'testfunctions' 0;'multipliers' 0;'current stepsize' 0};
    gds.numeric.LC = {'parameters' 1;'period' 1;'testfunctions' 0;'multipliers' 0;'current stepsize' 0;'PRC' 0;'dPRC' 0;'Input' 0};
    % XXXX
    gds.numeric.PD = {'parameters' 1;'period' 1;'testfunctions' 0;'multipliers' 0;'current stepsize' 0};
    % XXX
    gds.numeric.Hom = {'parameters' 1;'period' 1;'testfunctions' 0;'multipliers' 0;'current stepsize' 0};
    gds.numeric.HSN = {'parameters' 1;'period' 1;'testfunctions' 0;'multipliers' 0;'current stepsize' 0};
    gds.diagram = 'diagram';
    gds.parameters=[];gds.options.ActiveParams=[];
    gds.T=[];gds.eps0=[];gds.eps1=[];gds.extravec=[1 0 0]; gds.t=[]; gds.epsilon=[];
    gds.options.IgnoreSingularity=[];gds.plot2='';gds.plot3='';gds.PRC='';gds.dPRC=''; 
    gds.options.PRC = 0; gds.options.dPRC = 0;
         
%---------------------------------------------------------------------
function warning_box(stri1,stri2,handles)
global gds path_sys
button = questdlg('System already exist! Do you want to continue? If you press yes to continue, you will overwrite the existing system',...
'System already exist','Yes','No','No');
dir=path_sys;
if strcmp(button,'Yes')
   gds.system=stri1;
   delete(stri2);
elseif strcmp(button,'No')
   set(handles.name_system,'String','');
end

%---------------------------------------------------------------------------
function string=make_init
global gds;
string{1,1}='y0=[';
if (gds.dim>1)
    for i=1:(gds.dim-1)
        string{1,1}=strcat(string{1,1},'0,');
    end
end
string{1,1}=strcat(string{1,1},'0];');
string{2,1}=strcat('options = odeset(''Jacobian'',[],''JacobianP'',[],''Hessians'',[],''HessiansP'',[]);');   

%-----------------------------------------------------------------------------
function string=parsefile(file)
fid=fopen(file);
i=1;string=[];string{1,1}='empty';
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    string{i,1}='';
    while (~isempty(findstr(tline,'\'))||isempty(findstr(tline,']]')))
        tline=strrep(tline,'\','');
        if findstr(tline,'[[');
            tline=strrep(tline,'],','];');
            string{i,1}=tline;
        else
            tline=strrep(tline,'],','];');
            string{i,1}=strcat(string{i,1},tline);
        end
        tline=fgetl(fid);
    end
    tline=strrep(tline,'],','];');
    string{i,1}=strcat(string{i,1},tline);
    i=i+1;
end   
fclose(fid);


function line = symmat2line(symmat)
      [n,m] = size(symmat);
      line = '';
      for i = 1:n
	  row = '';
	  for j = 1:m
	      if (j ~= 1)
		  row = [row ' , ' char(symmat(i,j))];
	      else
		  row = char(symmat(i,j));
	      end
	  end
	  if (i ~= 1)
	      line = [line ' ; ' row];
	  else
	      line = row;
	  end
      end
      line  = ['[ ' line ' ]'];


% --- Executes during object creation, after setting all properties.
function name_system_CreateFcn(hObject, eventdata, handles)
% hObject    handle to name_system (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function coordinates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coordinates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function parameters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sys_Callback(hObject, eventdata, handles)
% hObject    handle to sys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sys as text
%        str2double(get(hObject,'String')) returns contents of sys as a double


% --- Executes during object creation, after setting all properties.
function sys_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- toGdsStruct('x,y,z')  ->      {'x'    [0], 'y'    [0], 'z'    [0]}
function result = toGdsStruct(str)
    items = strsplit(str, ',');
    result = [items', num2cell(zeros(length(items), 1))];

