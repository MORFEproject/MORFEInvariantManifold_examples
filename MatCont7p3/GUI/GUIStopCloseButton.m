classdef GUIStopCloseButton < handle

    
    properties
        handle
        endfunction
        eventlistener = [];
    end
    
    methods
        
        
        function obj = GUIStopCloseButton(parent , endfunction,  varargin )
        global sOutput
            obj.handle = uicontrol(parent , 'Style' , 'pushbutton'  , 'Callback' , @(src,ev) obj.buttonPushed(), varargin{:} ); 
            obj.eventlistener = sOutput.model.addlistener('stateChanged' , @(o,e) obj.updateLabel());
            set(obj.handle,'DeleteFcn' , @(o,e) obj.destructor());
            obj.endfunction = endfunction;
            obj.updateLabel();
        end
        
        
        function updateLabel(obj)
        global sOutput
            if (sOutput.currentstate == sOutput.ENDED)
              set(obj.handle , 'String' , 'Close');
            else
               set(obj.handle, 'String' , 'Stop');
           end
        end
        
        function buttonPushed(obj)
        global sOutput
           if (sOutput.currentstate == sOutput.ENDED)
               feval(obj.endfunction)
           else
               sOutput.performStop();
           end
        end
        
        
        function destructor(obj)
           %D disp(['Destructor: ' mfilename]);
           delete(obj.eventlistener);
           delete(obj);
        end
        
    end
    
end

