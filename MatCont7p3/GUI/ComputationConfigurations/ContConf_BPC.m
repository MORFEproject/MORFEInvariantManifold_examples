classdef ContConf_BPC < ContConf_LCCodim2
    properties
        pointlabel
        
    end
    
    methods
        function obj =  ContConf_BPC(pointlabel, initFunction)
            if nargin == 0
               pointlabel = 'BPC';
               initFunction = @init_BPC_BPC;
            end
            
            obj.curvedefinition = @branchpointcycle;
            obj.label = 'BPC';
            obj.defaultPointType = 'BPC';
            obj.pointlabel = pointlabel;
            obj.initFunction = initFunction;
            
            obj.testLabels = {};
            
        end
        function list = getGlobalVars(obj)
            list = {'cds', 'lds'};
            
        end
        function s = getLabel(obj)
            s = [getLabel@ContConf_LCCodim2(obj) '_' obj.label '_' obj.pointlabel];
        end
        
        
        function p = getPrioritynumber(~)
            p = 400;
        end
        
        function b = isAvailable(obj, settings)
            initialpoint = settings.IP;
            b = any(strcmp(initialpoint.getILabel(), obj.pointlabel));
        end
        function configureSettings(obj, settings)
            configureSettings@ContConf_LCCodim2(obj, settings);
            parammodel = settings.getSetting('parameters');
            parammodel.revive(settings, true, true); %active:true, branch:true
        end
        
        function msg = check(obj, settings)
            parametersmodel = settings.getSetting('parameters');
            activeParams = parametersmodel.getActive();
            bp = parametersmodel.getBranch();
            if (length(activeParams) ~= 3)
                msg = sprintf('You have to select exactly %s free parameters', num2str(3));
                return;
            end
            if (length(bp) ~= 1)
                msg = 'you have to select exactly one branch parameter';
                return;
            end
            msg = '';
            return;
        end    
        
    end
    
    methods(Access=protected)
        function [valid, msg] = sanityCheck(obj, handle, x0, param, activeParams, settings)
            msg = obj.check(settings);
            valid = isempty(msg);
        end
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
        
            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
  
            [x0,v0] = obj.initFunction(handle, IP.source.x,  sdata, activeParams,  settings.ntst, settings.ncol);
        end
    end
end
