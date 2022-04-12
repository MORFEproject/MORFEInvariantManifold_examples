classdef GUIMenuItem < handle

    
    properties
        handle
        model
        
        eventlistener = [];
        
        validator
        labelgetter
        callback
        
    end
    
    methods
        
        
        function obj = GUIMenuItem(parent,  callback, labelgetter, model, validator, eventdescr, varargin)
            
            if ~isa(labelgetter, 'function_handle')
               labelgetter = @() labelgetter; %if string is given, make a function that returns the string. 
            end
            
            obj.handle = uimenu(parent  , 'Callback' ,@(o, e) obj.executeCallback() , varargin{:});  
            
            obj.validator = validator;
            obj.labelgetter = labelgetter;
            obj.callback = callback;
            
            obj.updateLabel();
            set(obj.handle,'DeleteFcn' , @(o,e) obj.destructor());
            
            
            if ~isempty(model)
                obj.eventlistener = model.addlistener(eventdescr , @(srv,ev) obj.updateLabel() ); 
            end
            
        end
        
        function executeCallback(obj)
           if obj.validator()
              obj.callback(); 
           end
            
        end
        
        function updateLabel(obj)
            set(obj.handle,'Label', obj.labelgetter());
            set(obj.handle,'Enable' , CLbool2text(obj.validator()));
        end
        
        function destructor(obj)
           delete(obj.eventlistener);
           delete(obj);
        end
        
    end
    
end

