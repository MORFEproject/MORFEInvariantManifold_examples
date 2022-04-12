classdef CLSettingParameterSelect < CLSetting
    
    properties
        index
        parammodel
        tag
    end
    
    methods
        function obj = CLSettingParameterSelect(parammodel, index, tag, indexoffset)
            parameter = parammodel.parameters{index};
            obj = obj@CLSetting(['param_' parameter '_' tag], parammodel.get(tag, index), InputRestrictions.BOOL, 2, 2, indexoffset + index, '~~~~');
            obj.tag = tag;
            obj.index = index;
            obj.parammodel = parammodel;
            
            
            
        end
        
        function [valid, msg] = setValue(obj, newvalue)
            [valid, msg] = obj.validitycheck.validate(newvalue);
            if valid
                obj.parammodel.set(obj.tag, obj.index, newvalue);
            end
        end
        
        function value = getValue(obj)
            value = obj.parammodel.get(obj.tag, obj.index);
        end
        
        function newobj = copy(~, ~)
            newobj = [];
        end
        %{
       function b = isVisible(obj)
         b= obj.parammodel.isVisible(); 
       end
        %}
             function box = renderGUI(varargin)
          box = []; 
       end
    end
end
