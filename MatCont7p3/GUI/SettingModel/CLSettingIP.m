classdef CLSettingIP < handle
    properties
        currentpoint;
        currentpointdata;
        initiallist = {};
        initials = struct();
        
        source = [];
    end
    
    
    methods
        function obj = CLSettingIP(label)
            
            obj.initiallist = {...
                InitialType('Point','P'), ...
                '-', ...
                InitialType('Equilibrium','EP'), ...
                '-', ...
                InitialType('Limit cycle','LC'), ...
                '-', ...
                InitialType('Period Doubling','PD'), ...
                InitialType('Limit point of cycles','LPC'), ...
                InitialType('Neimark-Sacker','NS'), ...
                '-', ...
                InitialType('Hopf','H'), ...
                InitialType('Limit point','LP'), ...
                InitialType('Branch point','BP'), ...
                '-', ...
                InitialType('Cusp','CP'), ...
                InitialType('Bogdanov-Takens','BT'), ...
                InitialType('Generalized Hopf','GH'), ...
                InitialType('Zero-Hopf','ZH'), ...
                InitialType('Double Hopf','HH'), ...
                InitialType('Branch point of cycles','BPC'), ...
                '-', ...
                InitialType('Cusp Point of Cycle','CPC'), ...
                InitialType('Generalized Period Doubling','GPD'), ...
                '-', ...
                InitialType('Homoclinic to saddle','Hom'), ...
                InitialType('Homoclinic to saddle-node','HSN'), ...
                '-', ...
                InitialType('Non-central HomSN','NCH'), ...
                '-', ...
                InitialType('ConnectionSaddle','ConnecHom'), ...
                InitialType('ConnectionSaddleNode','ConnecHSN'), ...
                InitialType('HomotopySaddle','HTHom'), ...
                InitialType('HomotopySaddleNode','HTHSN'), ...
                '-', ...
                InitialType('ConnectionHet','ConnecHet'), ...
                InitialType('HomotopySaddleHet','HTHet'), ...
                InitialType('Heteroclinic','Het'), ...
                };
            
            for st = [obj.initiallist, {InitialType('<nothing>', '', 'EMPTY')}]
                st = st{1};
                if isobject(st)
                    obj.initials.(st.getILabel()) = st;
                end
            end
            
            if nargin > 0
                obj.setValue(label);
            end
            
        end
        
        function b = assignInternalLabel(obj, ilabel)
           obj.currentpoint.internallabel = ilabel;
           b = 1;
        end
        function b = assignLabel(obj, ilabel)
           obj.currentpoint.publiclabel = ilabel;
           b = 1;
        end        
        
        function b = setSource(obj, computation)
            obj.source = computation; b = 1;  %%% FIXME, CURVE Inflation problem
            
        end
      function value = getValue(obj)
            value = obj.currentpoint;
            
       end
       function [valid, msg] = setValue(obj, newvalue)
            valid = 1; msg = [];
            if isstruct(newvalue)
                label = strip(newvalue.label);
                if isfield(newvalue, 'internallabel') && ~isempty(newvalue.internallabel)
                   ilabel = strip(newvalue.internallabel); 
                else
                   ilabel = label;
                end
                    
               obj.currentpoint = InitialType(newvalue.msg, label, ilabel); 
               obj.currentpointdata = newvalue;
            else
                if isfield(obj.initials, newvalue)
                    obj.currentpoint = obj.initials.(newvalue);
                    obj.currentpointdata = [];
                else
                   valid = 0;
                   msg = 'Invalid Initial Point Type';
                end
            end
            
       end
  
       function forceValue(~, ~)
       end
       function id = getGroupID(~)
           id = 0;
       end
       function id = getSubGroupID(~)
           id = 0;
       end       
       function id = getItemID(~)
           id = 2;
       end      
              
       function row = getIDs(obj)
          row = [obj.getGroupID(), obj.getSubGroupID(), obj.getItemID()]; 
       end
       
       function s = toString(obj)
          s = sprintf('%s (%s)', obj.currentpoint.getName(), obj.currentpoint.getLabel());
       end
       function newobj = copy(obj, ~)
           newobj = CLSettingIP();
           newobj.currentpoint = obj.currentpoint;
           newobj.currentpointdata = obj.currentpointdata;
           newobj.source = obj.source;
       end
       function h = getHelpStr(~)
           h = 'no help';
       end
       function b = isVisible(~)
          b = 1; 
       end
       function setVisible(~, ~)
       end
       function b = isAdjusted(~)
           b = 1;
       end
       function box = renderGUI(varargin)
       box = [];
       end
       function t = getValueType(~)
            t = 'NONE';
       end
       function b = sanityCheck(obj, settings)
            b = 1;
       end
   end
    methods(Static)
       function install(settings, name, initval, restrict, catdata)
           if ~settings.activateParameter(name)
               settings.installParameter(name, CLSetting(name, initval, restrict, catdata(1), catdata(2), catdata(3), CLSettingHelp.getHelp(name) ));
           end
       end    
       
       function testPriority()
           cc = CompConfigurations();
           ip = CLSettingIP();
           settings = CLSettings();
           settings.addSetting('IP', CLSettingIP('EMPTY'));
           settings.addSetting('system', CLSystem('Systems/adapt2.mat'));
           for k = 1:length(ip.initiallist)
                if ~ischar(ip.initiallist{k})
                    point = ip.initiallist{k};
                    settings.IP.set(point.publiclabel);
                    options = strjoin( cellfun(@(x) sprintf('%s (%s)', x.toString(), num2str(x.getPrioritynumber())), cc.getConfigs(settings), 'un', 0) , ', ');
                    fprintf('%s: %s\n', point.publiclabel, options);
                end
           end
       end
       
    end
   
    
end
