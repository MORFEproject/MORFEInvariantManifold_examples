classdef CLSettingParameter < CLSetting
    
    properties
        index
        parammodel
    end
    
    methods
        function obj = CLSettingParameter(parammodel, index, restriction)
            if nargin < 3
               restriction =  InputRestrictions.NUM;
            end
            
            parameter = parammodel.parameters{index};
            obj = obj@CLSetting(['param_' parameter], parammodel.value(index), restriction, 2, 2, index, '~~~~');
            obj.index = index;
            obj.parammodel = parammodel;
            
        end
        
        function [valid, msg] = setValue(obj, newvalue)
            [valid, msg] = obj.validitycheck.validate(newvalue);
            if valid
                obj.parammodel.value(obj.index) = newvalue;
            end
        end
        
        function value = getValue(obj)
            value = obj.parammodel.value(obj.index);
        end
        
        function newobj = copy(~, ~)
            newobj = [];
        end
        function b = isVisible(obj)
            b= obj.parammodel.isVisible();
        end
        function box = renderGUI(varargin)
            box = [];
        end
    end
end
