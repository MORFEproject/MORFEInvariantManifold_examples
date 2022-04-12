classdef CLSettingCoordinate < CLSetting

    properties
       index
       coordmodel
    end
    
    methods
        function obj = CLSettingCoordinate(coordmodel, index, subcat)
            if nargin < 3
                subcat = 1;
            end
            coordinate = coordmodel.coordinates{index};
            obj = obj@CLSetting(coordinate, coordmodel.value(index), InputRestrictions.NUM, 2, subcat, index, '~~~~');
            obj.index = index;
            obj.coordmodel = coordmodel;
            
        end
        
        function [valid, msg] = setValue(obj, newvalue)
            [valid, msg] = obj.validitycheck.validate(newvalue);
            if valid
               obj.coordmodel.value(obj.index) = newvalue;  
            end
        end
        
       function value = getValue(obj)
           value = obj.coordmodel.value(obj.index);
       end        
      function newobj = copy(~, ~)
           newobj = [];
      end                
        
      function b = isVisible(obj)
         b= obj.coordmodel.isVisible(); 
      end
          
      
    end
end
