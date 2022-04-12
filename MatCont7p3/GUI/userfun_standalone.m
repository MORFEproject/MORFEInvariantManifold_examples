function varargout = userfun_standalone(varargin)
% USERFUN Application M-file for userfun_standalone.fig
%    FIG = USERFUN launch userfun_standalone GUI.
%    USERFUN('callback_name', ...) invoke the named callback.


% Last Modified by GUIDE v2.0 04-Sep-2002 11:20:31
global gds oldgds path_sys driver_window MC;
if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');
    oldgds = gds;
	% Generate a structure of handles to pass to callbacks, and store it. 
    gds.ok = false;
	handles = guihandles(fig);
	guidata(fig, handles);
    load_listbox1(handles);
	% Wait for callbacks to run and window to be dismissed:
	uiwait(fig);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);      
	end
end

% --------------------------------------------------------------------
function adbutton_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.adbutton.

global gds;
label = get(handles.label,'String');
name  = get(handles.name,'String');
if length(label)>2
    warndlg('A label exist of 2 characters')
    set(handles.label,'String',label(1:2));
elseif length(label)<2
    while length(label)<2
        label = sprintf('%s%c',label,char(32));
    end
    set(handles.label,'String',label(1:2));
end
if length(name)<1
    warndlg('You have to enter the name of the function');
    return;
end
if isfield(gds,'userfunction')%userfunctions already exists
    dim = size(gds.userfunction,2);
    namestr  = cellstr(char(gds.options.UserfunctionsInfo.name));
    labelstr = cellstr(char(gds.options.UserfunctionsInfo.label));
    gds.options.UserfunctionsInfo(dim+1).label = label;
    gds.options.UserfunctionsInfo(dim+1).name = name;
    gds.options.UserfunctionsInfo(dim+1).state = 1;
    gds.userfunction{dim+1} = get(handles.edituserfunction,'String');
    i = find(strcmp(namestr,name));
    if i
        button = questdlg('Userfunction already exist! Do you want to continue? If you press yes to continue, you will overwrite the existing userfunction',...
            'Userfunction already exist','Yes','No','No');
        if strcmp(button,'No')
            gds.userfunction(dim+1) = [];
            gds.options.UserfunctionsInfo(dim+1) = [];
            return;
        elseif strcmp(button,'Yes')
            gds.userfunction(dim+1)=[];
            gds.options.UserfunctionsInfo(dim+1) = [];
            gds.options.UserfunctionsInfo(i).label = label;
            gds.options.UserfunctionsInfo(i).name = name;
            gds.userfunction{i} = get(handles.edituserfunction,'String');
        end
    end
    if find(strcmp(labelstr,label))
        warndlg('There is another userfunction using this label: please choose another label!');
        set(handles.label,'String','');
        gds.userfunction(dim+1) = [];
        gds.options.UserfunctionsInfo(dim+1) = [];
    end
else
    gds.options.UserfunctionsInfo=[];
    gds.userfunction{1}=get(handles.edituserfunction,'String');
    gds.options.UserfunctionsInfo.label = label;
    gds.options.UserfunctionsInfo.name  = name;
    gds.options.UserfunctionsInfo.state = 1;
end
load_listbox1(handles);

  
% --------------------------------------------------------------------
function listbox1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.renamebutton.
global gds;
if isfield(gds,'userfunction')
    val=get(handles.listbox1,'Value');
    d=0;
    for k=1:size(gds.userfunction,2)
        if ~isempty(gds.userfunction{k})
            d=d+1;
            if d==val
                set(handles.name,'String',gds.options.UserfunctionsInfo(k).name);
                name_Callback(handles,[],handles);
                set(handles.label,'String',gds.options.UserfunctionsInfo(k).label);
                set(handles.edituserfunction,'String',gds.userfunction{k});
                break;
            end
        end
    end
else
    set(handles.name,'String','');
    name_Callback(handles,[],handles);
    set(handles.label,'String','');
    set(handles.edituserfunction,'String','');
end


% --------------------------------------------------------------------
function label_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.renamebutton.
global MC
label=get(handles.label,'String');
if length(label)>2
    warndlg('A label exist of 2 characters')
    set(handles.label,'String',label(1:2));
elseif length(label)<2
    while length(label)<2
        label=sprintf('%s%c',label,char(32));
    end
    set(handles.label,'String',label(1:2));
end


%list  = get(MC.mainwindow.initial_point,'children');
%tag   = get(list,'Tag');
%label = strcat(deblank(label),'_');
%i = strmatch(label,tag);
%if find(i)
%    warndlg('Labels should be unique and differ from labels of standard special points.')
%    set(handles.label,'String','');
%end

% --------------------------------------------------------------------
function name_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.renamebutton.
global gds
set(handles.adbutton,'Enable','on');
set(handles.updatebutton,'Enable','off');
name = get(handles.name,'String');
if isfield(gds,'userfunction')%userfunctions already exists
    if find(strcmp(cellstr(strvcat(gds.options.UserfunctionsInfo.name)),name))
        set(handles.updatebutton,'Enable','on');
        set(handles.adbutton,'Enable','off');
    end
end


% --------------------------------------------------------------------
function deletebutton_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.deletebuttton.
global gds
% set(handles.deletebutton,'Enable','off');

name  = get(handles.listbox1,'String');
val   = get(handles.listbox1,'Value');
label = get(handles.label,'String');
if isfield(gds,'userfunction')%userfunctions already exists
    dim = size(gds.userfunction,2);
    i = find(strcmp(cellstr(char(gds.options.UserfunctionsInfo.name)),name{val}));
    if i
        button = questdlg('Do you want to continue? If you press yes to continue, you will delete the existing userfunction',...
                'Delete userfunction ','Yes','No','No');          
        if strcmp(button,'Yes')
            if isfield(gds,'poincare_eq') && ~isempty(gds.poincare_eq)
                funct=func2str(gds.poincare_eq);
                if strcmp(funct,name{val})
                    gds.poincare_eq=[];
                    gds.poincare_do=0;
                    gds=rmfield(gds,'poincare_eq');
                end
            end    
            gds.userfunction{i} = '';
            gds.options.UserfunctionsInfo(i).state = 0;
        end
    end
end
load_listbox1(handles);

% --------------------------------------------------------------------
function updatebutton_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.updatebutton.
global gds;
nameu = get(handles.name,'String'); label = get(handles.label,'String');
if isfield(gds,'userfunction')%userfunctions already exists
    i = find(strcmp(cellstr(char(gds.options.UserfunctionsInfo.name)),nameu));
    gds.options.UserfunctionsInfo(i).label = label;
    gds.options.UserfunctionsInfo(i).name = nameu;
    gds.options.UserfunctionsInfo(i).state = 1;
    gds.userfunction{i}=get(handles.edituserfunction,'String');
end
load_listbox1(handles);

% --------------------------------------------------------------------
function cancelbutton_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.cancelbutton.
global gds oldgds path_sys;
gds  = oldgds;
file = fullfile(path_sys,gds.system);
save(file,'gds');
delete(handles.userfunfig);
if gds.open.figuur==1,starter;end
if gds.open.numeric_fig==1,numeric;end



%-----------------------------------------------------------------
function load_listbox1(handles)
global gds;
if isfield(gds,'userfunction') && ~isempty(char(gds.userfunction))
    str = [];d=0;
    for i=1:size(gds.userfunction,2)
        if ~isempty(gds.userfunction{i})
            d=d+1;
            str{d,1} = gds.options.UserfunctionsInfo(i).name;
            val=i;
        end
    end
    set(handles.listbox1,'String',cellstr(str),'Value',d);
    set(handles.name,'String',gds.options.UserfunctionsInfo(val).name);
    name_Callback(handles,[],handles);
    set(handles.label,'String',gds.options.UserfunctionsInfo(val).label);
    set(handles.edituserfunction,'String',gds.userfunction{val});
else
    set(handles.listbox1,'String','');
    set(handles.name,'String','');
    name_Callback(handles,[],handles);
    set(handles.label,'String','');
    set(handles.edituserfunction,'String','res=');

end


% --------------------------------------------------------------------
function okbutton_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.okbutton.
global gds path_sys MC;
name = get(handles.name,'String');
ad   = get(handles.adbutton,'Enable');
if ~isempty(name) && strcmp(ad,'on')
    adbutton_Callback(h,[],handles,[]);
end
updat = get(handles.updatebutton,'Enable');
if ~isempty(name) && strcmp (updat,'on')
    updatebutton_Callback(h,[],handles,[]);
end
if isfield(gds,'userfunction'),gds.options=contset(gds.options,'Userfunctions',1);
else gds.options=contset(gds.options,'Userfunctions',0);end
string_jac = '';string_jacp = '';string_hess = '';string_hessp = '';string_tensor3='';string_tensor4='';string_tensor5='';
dimp = size(gds.parameters,1);par='';pa='';
if ~isempty(dimp)
    par = cellstr((strcat(',',char(gds.parameters{:,1}))));
    par = strcat(par{:,1});
    pa  = par(2:end);
end
cor = '';
t = gds.time{1,1};
if (~isempty(gds.dim))
    cor = cellstr((strcat(',',char(gds.coordinates{:,1}))));
    cor = strcat(cor{:,1}); cor = cor(2:end);
end
if (~isempty(t))
    t=strcat(t,',');
else
    t='t,';
end
string_sys = cellstr(systems_standalone('replace_sys_input',gds.equations));
if strcmp(string_sys,'error')
    errordlg('The number of equations is not correct','Error');
    return
end
fwrite = strcat(gds.system,'.m');
fwrite = fullfile(path_sys,fwrite);
[fid_write,message] = fopen(fwrite,'w');
if fid_write==-1
    errordlg(message,'Error');
    return
end
fread = fullfile(path_sys,'standard.m');
[fid_read,message] = fopen(fread,'r');
if fid_read == -1
    errordlg(message,'Error');
    return
end
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
   
string_index = 3;
der_index = 1;
if sum(gds.der(2:end,der_index)) >= 1
    string_handles{string_index}=strrep(string_handles{string_index},'[]','@jacobian');
    string_handles{string_index+1}=strrep(string_handles{string_index+1},'[]','@jacobianp');
    string_index = string_index + 2;
end
der_index = 2;
if sum(gds.der(2:end,der_index)) >= 1
    string_handles{string_index}=strrep(string_handles{string_index},'[]','@hessians');
    string_handles{string_index+1}=strrep(string_handles{string_index+1},'[]','@hessiansp');
    string_index = string_index + 2;    
end
for der_index = 3:5
    if sum(gds.der(2:end,der_index)) >= 1
            string_handles{string_index} = strrep(string_handles{string_index},'[]', sprintf('@der%i', der_index));    
    end
    string_index = string_index + 1;   
end
            

                
if ~isempty(gds.options.UserfunctionsInfo)
    siz = size(gds.options.UserfunctionsInfo,2);
    for i = 1:siz
        string_handles{9+i,1}= sprintf('out{%d}= @%s;',9+i,gds.options.UserfunctionsInfo(i).name);
    end
else siz=0;end



fprintf(fid_write, strrep(fgetl(fid_read),'odefile',gds.system));
fprintf(fid_write, '\n');

for i=1:9+siz
    fprintf(fid_write,'%s\n',string_handles{i,1});
end
fprintf(fid_write, '%s', gds.filecontent);


if ~isempty(gds.options.UserfunctionsInfo)    
   for i=1:size(gds.options.UserfunctionsInfo,2)
       res=0;
       if isfield(gds,'userfunction') && ~isempty(gds.userfunction{i})
           str_user = systems_standalone('replace_token',cellstr(gds.userfunction{i}));
       else 
           str_user=cellstr('res=');
       end
       [userline, ~, newpar] = renameforsym('', '', strip(par, ','));
       [~,~,newpar] = renameforsym('', '', strip(par, ','));
       newpar = strcat(',', newpar);
       hs1 = sprintf('function userfun%d=%s(t,kmrgd%s)',i,gds.options.UserfunctionsInfo(i).name,newpar);
       fprintf(fid_write,'%s\n',hs1);
       hs1 = sprintf('userfun%d',i);
       dim = size(str_user,1);
       for j = 1:dim
           userline = renameforsym(str_user{j}, '', strip(par, ','));
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
fclose(fid_read);
fclose(fid_write);
file = fullfile(path_sys,gds.system);
delete(handles.userfunfig);
save(file,'gds');
gds.ok = true;
rehash;
