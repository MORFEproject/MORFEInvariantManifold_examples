classdef GUIEvalPlotOpsBox < handle

    properties
        handle
        setter
        getter

        eventlistener = [];
        validator;
    end
    
    methods
        
        
        function obj = GUIEvalPlotOpsBox(parent , setter , getter , model, event , varargin )
            
            celllist = getter();
            obj.handle = uicontrol(parent , 'Style' , 'edit' , 'Unit' , 'Pixels'  , 'String' , obj.celllist2string(celllist) , 'Callback' , @(src,ev) obj.newValue()  , varargin{:} , 'TooltipString', 'enter MATLAB plot options');  
            set(obj.handle,'DeleteFcn' , @(o,e) obj.destructor());
            obj.setBackgroundOK();
            obj.setter = setter;
            obj.getter = getter;

            
            if ~isempty(model)
                obj.eventlistener = model.addlistener(event , @(srv,ev) obj.settingChanged()); 
            end
            
            obj.validator = @(x) GUIEvalPlotOpsBox.testIfValidPlotOptions(x);
            
        end
        function setTextMode(obj)
           obj.validator = @(x) GUIEvalPlotOpsBox.testIfValidTextOptions(x);
        end
        
        
        function newValue(obj)
             string = get(obj.handle ,'String');
             if (isempty(strtrim(string)))   
                obj.settingChanged(); %%used to restore default value; 
                return;
             else
            
                try 
                     x = evalin('base' , ['{' string  '}' ]);
                     if (iscell(x) && (min(size(x)) == 1))
                        [bool , errormsg] = obj.validator(x);
                        if (bool)
                            obj.setBackgroundOK();
                            obj.setter(x);
                        else
                            fprintf('Error while evaluating entry: PLOT option error, ');
                            obj.setBackgroundERROR();
                            disp(errormsg);
                            obj.informError('Error while evaluating entry:', errormsg);
                        end
                     else
                         obj.setBackgroundERROR();
                         fprintf('Error while evaluating entry: ');
                         fprintf('Entry is not a Cell list\n'); 
                         obj.informError('Error while evaluating entry');
                         
                         
                     end
                catch error
                     obj.setBackgroundERROR();
                     fprintf('Error while evaluating entry: ');
                     disp(error.message);
                     obj.informError('Error while evaluating entry:', error.message);
                end
            end
        end
        
        function settingChanged(obj)
           set(obj.handle , 'String' , obj.celllist2string(obj.getter())); 
           obj.setBackgroundOK();
        end
        
        function destructor(obj)
            delete (obj.eventlistener);
            delete(obj);
        end

        function setBackgroundOK(obj)
           set(obj.handle, 'BackgroundColor' , [0.9 1 0.9]);
        end

        function setBackgroundERROR(obj)
           set(obj.handle, 'BackgroundColor' ,  [1    0.3    0.3]);
        end
        function e = Extent(obj)
           e = obj.handle.Extent; 
        end
        function set(obj, varargin)
           set(obj.handle, varargin{:}); 
        end                      
    end
    
    methods(Static)
        function informError(varargin)
            string = varargin{1};
            for i = 2:length(varargin)
               string = [string , newline , varargin{i}]; 
            end
            errordlg(string , 'Syntax error');
            
        end
        
        
        function s = celllist2string(cellList)
            s = strjoin(cellfun(@(x) mat2str(x), cellList, 'un', 0), ', ');
        end
        
        function [bool , errormsg] = testIfValidPlotOptions(plotops)
            try
                fig = figure('Visible' , 'off');
                ax = axes('Parent' ,fig , 'Visible' , 'off');
                line([0 1] , [0 1] , [0 1] , 'Parent' , ax , plotops{:});
                delete(ax);
                delete(fig);
                bool = true;
                errormsg = '';
            catch error
                errormsg = error.message;
                bool = false;
            end
        end
        function [bool , errormsg] = testIfValidTextOptions(plotops)
         
            try
                fig = figure('Visible' , 'off');
                ax = axes('Parent' ,fig , 'Visible' , 'off');
                text([0 1] , [0 1] , [0 1] , 'dummy',  'Parent' , ax , plotops{:});
                delete(ax);
                delete(fig);
                bool = true;
                errormsg = '';
            catch error
                errormsg = error.message;
                bool = false;
            end
        end        
    end
end

