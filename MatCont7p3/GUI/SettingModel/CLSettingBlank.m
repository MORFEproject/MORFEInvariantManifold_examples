classdef CLSettingBlank < CLSetting
    properties
        blankmessage
    end
    
    methods
        function obj = CLSettingBlank(displayname, initvalue, checker, groupid, subgroupid, itemid, help, blankmessage)
            obj = obj@CLSetting(displayname, initvalue, checker, groupid, subgroupid, itemid, help);
            obj.blankmessage = blankmessage;
        end
        
        function [valid, msg] = setValue(obj, newvalue)
            if isempty(newvalue)
                obj.value = [];
                valid = true; msg = '';
            else
                [valid, msg] = setValue@CLSetting(obj, newvalue);
            end
        end
        
        function s = toString(obj)
            if isempty(obj.value)
                s = obj.blankmessage;
            else
                s = toString@CLSetting(obj);
            end
        end
        function newobj = copy(obj, newsettings)
            newobj = CLSettingBlank(obj.displayname, obj.value, obj.validitycheck, obj.groupid, obj.subgroupid, obj.itemid, obj.help, obj.blankmessage);
            obj.copyOver(newobj);
        end
        
        
    end
    
end