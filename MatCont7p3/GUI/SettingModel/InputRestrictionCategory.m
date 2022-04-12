classdef InputRestrictionCategory

    properties
       categories
    end
    methods
        function obj =  InputRestrictionCategory(categories)
           obj.categories = categories;
        end
        
        function [passed, msg] = validate(obj, value)
            passed = ~isempty(value) && any(contains(obj.categories, value));
            if passed
               msg = '';  
            else
               msg = ['one of these values should be set: ', sprintf('%s, ', obj.categories{:})];
               msg = msg(1:end-2);
            end
        end
        
        function value = adjust(~, value)
        end
        
        function str = toString(~, value)
            str = num2str(value);
        end
 
        function t = getType(~)
            t = 'CATEGORY';
        end

        function r = getGUIRenderer(~)
                r = @GUISettingPopupmenu;
        end
    end

end
