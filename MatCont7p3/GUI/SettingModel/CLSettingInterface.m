classdef CLSettingInterface < handle
    
    
    methods
        % % % CLSetting INTERFACE % % %
        function value = getValue(obj); value = obj; end
        function [valid, msg] = setValue(~, ~); valid = 0; msg = 'not allowed'; end
        function id = getGroupID(~); id = 0; end
        function id = getSubGroupID(~); id = 0; end
        function id = getItemID(~); id = 0; end
        function s = toString(~); s = '<settingsinterface>'; end
        function h = getHelpStr(~); h = ''; end
        function b = isVisible(~); b = 0; end
        function setVisible(~, ~); end
        function row = getIDs(obj); row = [obj.getGroupID(), obj.getSubGroupID(), obj.getItemID()]; end
        function t = getValueType(~); t = 'NONE'; end
        function b = sanityCheck(~, ~); b = 1; end
        
        function newobj = copy(~, ~)
            assert(0, 'copy needs to be implemented\n'); newobj = [];
        end
        
    end
end