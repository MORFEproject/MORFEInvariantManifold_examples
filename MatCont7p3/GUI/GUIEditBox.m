classdef GUIEditBox < handle
    
    properties
        handle
        OKColor;
        getter;
        setter;
        validator;
    end
    
    methods
        
         function obj = GUIEditBox(parent, getter, setter, validator, varargin )
             

            obj.getter = getter;
            obj.setter = setter;
            obj.validator = validator;
             
            obj.handle = uicontrol(parent , 'Style' , 'edit'  , 'Unit' , 'Pixels'  , 'String' , '' , 'Callback' , @(src,ev) obj.newValue()  , varargin{:} );  
            obj.OKColor = get(obj.handle, 'BackgroundColor');
            set(obj.handle,'DeleteFcn' , @(o,e) obj.destructor());
            obj.setBackgroundOK();
            obj.syncval();
         end
        
         function syncval(obj)
             set(obj.handle, 'String', num2str(obj.getter() , '%.8g'))
             
         end
         
        function newValue(obj)
            string =  get(obj.handle ,'String');

            try
                if ~isempty(string)
                    x = evalin('base' , string);
                else
                    x = []; 
                end
                [valid, errormsg] = obj.validator.validate(x);
                
                if (~valid)
                    obj.performErrorDisplay();
                    fprintf(2, sprintf('[%s] ERROR(%s): %s, value: %s\n\n',datetime('now', 'format', 'HH:mm:ss'), '', errormsg, string));
                    obj.syncval();
                else
                    obj.setter(x);
                end
                
            catch error
          
                fprintf(2, sprintf('[%s] ERROR(%s): %s: %s, value: %s\n\n',datetime('now', 'format', 'HH:mm:ss'), '', error.message, string));
                obj.performErrorDisplay();
                obj.syncval();
                
            end
       
  
        end
        function performErrorDisplay(obj)
            obj.setBackgroundERROR();
            pause(0.3)
            obj.setBackgroundOK();
            
        end
        function setBackgroundOK(obj)
           set(obj.handle, 'BackgroundColor' , [1 1 1]);
           set(obj.handle, 'BackgroundColor' , obj.OKColor);
        end

        function setBackgroundERROR(obj)
           set(obj.handle, 'BackgroundColor' ,  [1    0.3    0.3]);
        end
                       
        
        function destructor(obj)
            delete(obj);
        end
        
        function e = Extent(obj)
           e = obj.handle.Extent; 
        end
        function set(obj, varargin)
           %disp(varargin{2});
           set(obj.handle, varargin{:}); 
        end
    end
    
end

