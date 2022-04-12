function varargout= set_option(varargin)
% SET_OPTION Application M-file 
% last change on 13/11/02 11:38  
global gds path_sys MC;
if nargin ==1
    h=varargin{1};
else
    h = gcbo;
end
tag  = get(h,'Tag');
file = fullfile(path_sys,gds.system);
switch tag
case 'MinStepsize'     
      val = str2num(get(h,'String'));      
      if isnumeric(val)         
          gds.options = contset(gds.options, 'MinStepsize' ,val);
          save(file,'gds');
      end
case 'MaxStepsize'    
    val = str2num(get(h,'String'));
    if isnumeric(val)
        gds.options = contset(gds.options, 'MaxStepsize' ,val);
        save(file,'gds');
    end
case  'InitStepsize'    
    val = str2num(get(h,'String'));
    if isnumeric(val)
        gds.options = contset(gds.options, 'InitStepsize' ,val);    
        save(file,'gds');
    end
case 'FunTolerance' 
    val = str2num(get(h,'String'));
    if isnumeric(val)
        gds.options = contset(gds.options, 'FunTolerance' ,val);
        save(file,'gds');
    end
case 'VarTolerance'    
    val = str2num(get(h,'String'));
    if isnumeric(val)
        gds.options = contset(gds.options, 'VarTolerance',val);
        save(file,'gds');
    end
case 'TestTolerance'   
    val = str2num(get(h,'String'));
    if isnumeric(val)
        gds.options = contset(gds.options, 'TestTolerance' ,val);
        save(file,'gds');
    end
case 'MaxNewtonIters'    
    val = abs(round(str2num(get(h,'String'))));
    if val >= 0
        set(h,'String',val);
        gds.options = contset(gds.options, 'MaxNewtonIters' ,val);
        save(file,'gds');
    else
        set(h,'String',3);
    end
case 'MaxCorrIters'    
    val = abs(round(str2num(get(h,'String'))));
    if val>=0
        set(h,'String',val);
        gds.options = contset(gds.options, 'MaxCorrIters' ,val);
        save(file,'gds');
    else
        set(h,'String',10);
    end
case 'MaxTestIters'   
    val = abs(round(str2num(get(h,'String'))));
    if val>=0
        set(h,'String',val);
        gds.options = contset(gds.options, 'MaxTestIters' ,val);
        save(file,'gds');
    else
        set(h,'String',10);
    end
case 'MaxNumPoints'    
    val = abs(round(str2num(get(h,'String'))));
    if val>0
        set(h,'String',val);
        gds.options = contset(gds.options, 'MaxNumPoints' ,val);
        save(file,'gds');
    else
        set(h,'String',300);
    end
case 'CheckClosed'   
    val = abs(round(str2num(get(h,'String'))));
    if val>=0
        set(h,'String',val);
        gds.options = contset(gds.options, 'CheckClosed' ,val);
        save(file,'gds');
    else
        set(h,'String',10);
    end
case 'Increment'   
    val = str2num(get(h,'String'));
    if isnumeric(val)
        gds.options = contset(gds.options, 'Increment' ,val);
        save(file,'gds');
    end
case 'amplitude'   
    val = str2num(get(h,'String'));
    if isnumeric(val)
        gds.amplitude = val;
        save(file,'gds');
    end
case 'Adapt'   
    val=abs(round(str2num(get(h,'String'))));
    if  val >= 0
        set(h,'String',val);
        gds.options=contset(gds.options, 'Adapt' ,val);
        save(file,'gds');
    else
        set(h,'String',1);
    end
case 'ncol'
    val=abs(round(str2num(get(h,'String'))));
    if val >1 & val < 8
        set(h,'String',val);
        gds.discretization.ncol = val;
        save(file,'gds');
    else
        set(h,'String',4);
    end
case 'ntst'
    val = abs(round(str2num(get(h,'String'))));
    if val > 0
        set(h,'String',val);
        gds.discretization.ntst = val;
        save(file,'gds');
    else
        set(h,'String',20);
    end
case 'Multipliers'
    popup = get(h,'Value');    
    switch popup
    case 1
        val = 1;       
        gds.options = contset(gds.options, 'Multipliers' ,val);
        save(file,'gds');                        
    case 2
        val = 0;
        gds.options = contset(gds.options, 'Multipliers' ,val);   
        save(file,'gds');
    end    
case 'Eigenvalues'
    popup = get(h,'Value');    
    switch popup
    case 1
        val = 1;       
        gds.options = contset(gds.options, 'Eigenvalues' ,val); 
        save(file,'gds');                        
    case 2
        val = 0;
        gds.options = contset(gds.options, 'Eigenvalues' ,val);  
        save(file,'gds');
    end    

case 'ActiveParams'
    m = get(h,'Value');
    user = get(h,'UserData');
    num=user.num;
    if isempty(gds.options.ActiveParams)
        str = '[]';
        st  = [];
    else
        str = mat2str(gds.options.ActiveParams);
        st  = gds.options.ActiveParams;
    end
    switch m
    case 0
        x = find(st==str2num(num));
        if isempty(x)
            return
        else
            st(x) = '';str = mat2str(st);
            val = sort(str2num(str2mat(str)));
            gds.options = contset(gds.options, 'ActiveParams' ,val);
            save(file,'gds');
        end            
    case 1
        user = get(h,'UserData');num=user.num;
        if isempty(gds.options.ActiveParams)
            str = '[]';
        else
            str = mat2str(gds.options.ActiveParams);
        end
        if (size(gds.options.ActiveParams,2)==1)
           val = sort(str2num(str2mat(sprintf('[%s%c%s]',str,char(32),num))));
           gds.options = contset(gds.options, 'ActiveParams' ,val);   
           save(file,'gds');
        else
            str(end) = '';
            val = sort(str2num(str2mat(sprintf('%s%c%s]',str,char(32),num))));
            gds.options = contset(gds.options, 'ActiveParams' ,val);
            save(file,'gds');
        end
    end
case 'ActiveUParams'
    m = get(h,'Value');
    user = get(h,'UserData');
    num=user.num;
    if isempty(gds.options.ActiveUParams)
        str = '[]';
        st  = [];
    else
        str = mat2str(gds.options.ActiveUParams);
        st  = gds.options.ActiveUParams;
    end
    switch m
    case 0
        x = find(st==str2num(num));
        if isempty(x)
            return
        else
            st(x) = '';str = mat2str(st);
            val = sort(str2num(str2mat(str)));
            gds.options = contset(gds.options, 'ActiveUParams' ,val);
            save(file,'gds');
        end            
    case 1
        user = get(h,'UserData');num=user.num;
        if isempty(gds.options.ActiveUParams)
            str = '[]';
        else
            str = mat2str(gds.options.ActiveUParams);
        end
        if (size(gds.options.ActiveUParams,2)==1)
           val = sort(str2num(str2mat(sprintf('[%s%c%s]',str,char(32),num))));
           gds.options = contset(gds.options, 'ActiveUParams' ,val);   
           save(file,'gds');
        else
            str(end) = '';
            val = sort(str2num(str2mat(sprintf('%s%c%s]',str,char(32),num))));
            gds.options = contset(gds.options, 'ActiveUParams' ,val);
            save(file,'gds');
        end
    end    
    
case 'ActiveSParams'
    m = get(h,'Value');
    user = get(h,'UserData');
    num=user.num;
    if isempty(gds.options.ActiveSParams)
        str = '[]';
        st  = [];
    else
        str = mat2str(gds.options.ActiveSParams);
        st  = gds.options.ActiveSParams;
    end
    switch m
    case 0
        x = find(st==str2num(num));
        if isempty(x)
            return
        else
            st(x) = '';str = mat2str(st);
            val = sort(str2num(str2mat(str)));
            gds.options = contset(gds.options, 'ActiveSParams' ,val);
            save(file,'gds');
        end            
    case 1
        user = get(h,'UserData');num=user.num;
        if isempty(gds.options.ActiveSParams)
            str = '[]';
        else
            str = mat2str(gds.options.ActiveSParams);
        end
        if (size(gds.options.ActiveSParams,2)==1)
           val = sort(str2num(str2mat(sprintf('[%s%c%s]',str,char(32),num))));
           gds.options = contset(gds.options, 'ActiveSParams' ,val);   
           save(file,'gds');
        else
            str(end) = '';
            val = sort(str2num(str2mat(sprintf('%s%c%s]',str,char(32),num))));
            gds.options = contset(gds.options, 'ActiveSParams' ,val);
            save(file,'gds');
        end
    end    
    
case 'ActiveSParam'
    m = get(h,'Value');
    user = get(h,'UserData');
    num=user.num;
    if isempty(gds.options.ActiveSParam)
        str = '[]';
        st  = [];
    else
        str = mat2str(gds.options.ActiveSParam);
        st  = gds.options.ActiveSParam;
    end
    switch m
    case 0
        x = find(st==str2num(num));
        if isempty(x)
            return
        else
            st(x) = '';str = mat2str(st);
            val = sort(str2num(str2mat(str)));
            gds.options = contset(gds.options, 'ActiveSParam' ,val);
            save(file,'gds');
        end            
    case 1
        user = get(h,'UserData');num=user.num;
        if isempty(gds.options.ActiveSParam)
            str = '[]';
        else
            str = mat2str(gds.options.ActiveSParam);
        end
        if (size(gds.options.ActiveSParam,2)==1)
           val = sort(str2num(str2mat(sprintf('[%s%c%s]',str,char(32),num))));
           gds.options = contset(gds.options, 'ActiveSParam' ,val);   
           save(file,'gds');
        else
            str(end) = '';
            val = sort(str2num(str2mat(sprintf('%s%c%s]',str,char(32),num))));
            gds.options = contset(gds.options, 'ActiveSParam' ,val);
            save(file,'gds');
        end
    end        
    
    
case 'IgnoreSingularity' 
    
    popup = get(h,'Value');
    user  = get(h,'UserData');
    num=user.num;
    if isempty(gds.options.IgnoreSingularity)
        str='[]';
        st = [];
    else
        str = mat2str(gds.options.IgnoreSingularity);
        st  = gds.options.IgnoreSingularity;
    end
    switch popup
    case 1 %Yes
        x = find(st==str2num(num));
        if isempty(x)
            return
        else
            st(x) = [];str=mat2str(st);
            val = str2num(str) ;      
            gds.options = contset(gds.options, 'IgnoreSingularity' ,val);
            save(file,'gds');
        end            
    case 2 %no
       x = find(st==str2num(num)); 
       if ~isempty(x)
            return
        end
        if (size(gds.options.IgnoreSingularity,2)==1)
           val = str2num(str2mat(sprintf('[%s%c%s]',str,char(32),num)));
           gds.options = contset(gds.options, 'IgnoreSingularity' ,val);   
           save(file,'gds');
        else
            str(end) = '';
            val = str2num(str2mat(sprintf('%s%c%s]',str,char(32),num)));
            gds.options = contset(gds.options, 'IgnoreSingularity' ,val);   
            save(file,'gds');
        end
    end
    if all(size(gds.options.IgnoreSingularity)==[1 0])
        gds.options=contset(gds.options,'IgnoreSingularity',[]);
    end
    feval(gds.gui_singularities);
case 'UserfunctionsInfo' 
    popup = get(h,'Value');
    user  = get(h,'UserData');num=str2num(user.num);
    switch popup
    case 1 %Yes
        gds.options.UserfunctionsInfo(num).state=1;  
        gds.options = contset(gds.options,'Userfunctions',1);
    case 2 %no
        gds.options.UserfunctionsInfo(num).state=0;
        for i=1:size(gds.options.UserfunctionsInfo,2)
            if (gds.options.UserfunctionsInfo(i).state==1)
                return;
            end
        end
        gds.options = contset(gds.options,'Userfunctions',0);
    end
case 'DOfunction'
    popup = get(h,'Value');
    switch popup
        case 1 %none
            gds.poincare_do=0;gds.poincare_eq='';
        case 2 % add function
             userfun;
        otherwise
            d=0;
            for k=1:size(gds.userfunction,2)
                if ~isempty(gds.userfunction{k})
                    d=d+1;
                    if d==(popup-2)
                        gds.poincare_do=1;
                        funhandle = feval(gds.system);
                        gds.poincare_eq=funhandle{9+k};
                        gds.integrator.Crossection=gds.options.UserfunctionsInfo(k).name;
                        break;
                    end
                end
            end
     end     
 case 'nr_DOfunction'
    popup = get(h,'Value');
    switch popup 
        case 1% all
            gds.integrator.nr_crossection = 0;
        case 2%+1 if only zeros where the event function is increasing
            gds.integrator.nr_crossection = 1;
        case 3
            gds.integrator.nr_crossection = -1;
    end

% XXXX
case 'PRC'
    popup = get(h,'Value');    
    switch popup
    case 1 
        if length(MC.PRC) == 0
            MC.PRC =[];
        end
        if ~isempty(MC.PRC)
            close(MC.PRC);
        end  
        gds.options.PRC = 0;
        MC.PRC = [];                
    case 2
        val = 1;
        gds.options = contset(gds.options, 'PRC' ,val);  
        PRC_plot;
    end    

case 'dPRC'
    popup = get(h,'Value');    
    switch popup
    case 1
        if length(MC.dPRC) == 0
            MC.dPRC =[];
        end   
        if ~isempty(MC.dPRC)
            close(MC.dPRC);
        end  
        gds.options.dPRC = 0;
        MC.dPRC = [];          
    case 2
        val = 1;
        gds.options = contset(gds.options, 'dPRC' ,val);  
        dPRC_plot;
    end    
    
case 'Input'
    val = abs(str2num(get(h,'String')));
    if ~isempty(val)
        % double value
        if val == 0
            set(h,'String',0);
            gds.options = contset(gds.options,'Input',val);
        else
            gds.options = contset(gds.options,'Input',val);
            save(file,'gds');
        end
    else
        % function handle
        val = str2func(get(h,'String'));
        gds.options = contset(gds.options,'Input',val);
        save(file,'gds');
    end
    
case 'Tpar'
    m = get(h,'Value');
    gds.extravec(1) = m;
    
case 'eps0par'
    m = get(h,'Value');
    gds.extravec(2) = m;
    
case 'eps1par'
    m = get(h,'Value');
    gds.extravec(3) = m;
    
case 'T'
    val = str2num(get(h,'String'));
    if isnumeric(val)
        gds.T = val;
        save(file,'gds');
    end
    
case 'TToleranceedit'
    val = str2num(get(h,'String'));
    if isnumeric(val)
        gds.TTolerance = val;
        save(file,'gds');
    end      

case 'amplitudeedit'
    val = str2num(get(h,'String'));
    if isnumeric(val)
        gds.amplitude = val;
        save(file,'gds');
   end         
    
case 'epsedit'
    val = str2num(get(h,'String'));
    if isnumeric(val)
        gds.eps = val;
        save(file,'gds');
   end     

case 'eps1tol'
    val = str2num(get(h,'String'));
    if isnumeric(val)
        gds.eps1tol = val;
        save(file,'gds');
    end
   
% XXXX

otherwise
    errordlg(['Continuer option Error. Unknown tag: ' tag]);
end
