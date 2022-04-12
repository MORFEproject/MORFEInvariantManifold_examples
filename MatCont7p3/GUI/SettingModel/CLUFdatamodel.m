classdef CLUFdatamodel < handle
    properties
        ufdata = [];
    end
    
    methods
        function obj = CLUFdatamodel(ufdata)
            obj.ufdata = ufdata;
        end
        
        function value = getValue(obj)
            value = obj;
        end
        
        function [valid, msg] = setValue(obj, newvalue)
            valid = 0;
            msg = 'not allowed';
        end
        
        function id = getGroupID(obj)
            id = 0;
        end
        function id = getSubGroupID(obj)
            id = 0;
        end
        function id = getItemID(obj)
            id = 1;
        end
        
        function s = toString(obj)
            s = '<ufdata>';
        end
        
        function newobj = copy(obj, newsetting)
            newobj = CLUFdatamodel(obj.ufdata);
        end
        function h = getHelpStr(obj)
            h = 'no help';
        end
        function b = isVisible(obj)
            b = 1;
        end
        function setVisible(obj, bool)
        end
        function row = getIDs(obj)
            row = [obj.getGroupID(), obj.getSubGroupID(), obj.getItemID()];
        end
        function box = renderGUI(varargin)
            box = [];
        end
       function t = getValueType(~)
            t = 'NONE';
       end
       
       
       function b = sanityCheck(obj, settings)
            b = obj.internalSanityCheck(settings);
           
            if ~b
                uflabels = {obj.ufdata.label};
                for i = 1:length(uflabels)
                    uf = uflabels{i};
                    settings.removeSetting(['uf_', uf]);
                end
            end
            
       end
       
       function b = internalSanityCheck(obj, settings)
           system = settings.system;
           if isempty(system); b = 0; return; end
           
           
           sysuflabels = {system.ufdata.label};
           uflabels = {obj.ufdata.label};
          
          if length(sysuflabels) ~= length(uflabels); b = 0; return; end
          
          for k = 1:length(uflabels)
                if ~strcmp(sysuflabels{k}, uflabels{k})
                    b = 0;
                    return
                end
          end
          b = 1;
          
       end           
       
    end
    
end
