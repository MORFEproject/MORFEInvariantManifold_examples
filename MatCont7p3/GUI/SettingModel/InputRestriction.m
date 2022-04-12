classdef InputRestriction

    properties
       validator
       adjuster
       tostring
       failmsg
       guirenderer
       type
    end
 
       
    methods
        function obj =  InputRestriction(validator, adjuster, tostring, failmsg, guirenderer, type)
           obj.validator = validator;
           obj.adjuster = adjuster;
           obj.tostring = tostring;
           obj.failmsg = failmsg;
            
           if nargin < 5; obj.guirenderer =  @GUISettingEditBox; else; obj.guirenderer = guirenderer; end
           if nargin < 6; obj.type = 'VALUE'; else; obj.type = type; end
        end
        function r = getGUIRenderer(obj)
                r = obj.guirenderer;
        end
        function [passed, msg] = validate(obj, value)
            passed = obj.validator(value);
            if passed
               msg = '';  
            else
               msg = obj.failmsg;
            end
            
        end
        
        function value = adjust(obj, value)
           value = obj.adjuster(value); 
        end
        
        function str = toString(obj, value)
            str = obj.tostring(value);
        end
        function t = getType(obj)
            t = obj.type;
        end
    end
    
    methods(Static)
        function s = vector2string(x)
            s = ['[ ' sprintf('%.15g, ',x) ']'];
            s(length(s) - 2) = []; %%removes extra ',' at the end:  example  [ 5, 0.234324, 1e-12, ] --> [ 5, 0.234324, 1e-12 ]
        end
        function str = bool2string(b)
          if b
              str = 'true';
          else
              str = 'false';
          end

        end

    end

end
