classdef CLSettingUserFunction < CLSetting
    properties
       UFdatamodel 
       UFindex
    end
    
    methods
        function obj = CLSettingUserFunction(UFdatamodel, UFindex)
            
          
            name = sprintf('%s (%s)', UFdatamodel.ufdata(UFindex).label, UFdatamodel.ufdata(UFindex).name);
            state = UFdatamodel.ufdata(UFindex).state;
            obj = obj@CLSetting(name, state, InputRestrictions.BOOL, 2, 7, UFindex, '~~~');
            obj.UFdatamodel = UFdatamodel;
            obj.UFindex = UFindex;         
            
        end
        
        function value = getValue(obj)
           value = obj.UFdatamodel.ufdata(obj.UFindex).state;
           obj.value = value;
        end
        function [valid, msg] = setValue(obj, newvalue)
           [valid, msg] = obj.validitycheck.validate(newvalue);
           if valid
              obj.UFdatamodel.ufdata(obj.UFindex).state = newvalue; 
              obj.value = newvalue;
           end            
        end
        
        function newobj = copy(obj, newsettings)
            newUFdatamodel = newsettings.getSetting('userfunctions');
            if ~isempty(newUFdatamodel)
                newobj = CLSettingUserFunction(newUFdatamodel, obj.UFindex);
                obj.copyOver(newobj);
            else
                newobj = [];
                %CLdebugprint('Userfunction entry has died\n');
            end
            
        end
        
        
    end
end
