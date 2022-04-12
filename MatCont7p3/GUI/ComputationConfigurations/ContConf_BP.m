classdef ContConf_BP < ContConf_EP
    properties
        pointlabel
        
    end
    
    methods
        function obj =  ContConf_BP(pointlabel, initFunction)
            if nargin == 0
               pointlabel = 'BP';
               initFunction = @init_BP_BP;
            end
            
            
            obj.curvedefinition = @branchpoint;
            obj.label = 'BP';
            obj.nrActive = 3; 
            obj.defaultPointType = 'BP';
            obj.pointlabel = pointlabel;
            obj.initFunction = initFunction;
            
            obj.testLabels = {};
            
            
        end
        function list = getGlobalVars(obj)
            list = {'cds', 'bpds'};
            
        end
        
        function configureSettings(obj, settings)
            configureSettings@ContConf_EP(obj, settings);

            paramsetting = settings.getSetting('parameters');
            paramsetting.revive(settings, true, true); %active:true, branch:true    
            
      
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf_EP(obj) '_' obj.pointlabel];
        end
        
        function p = getPrioritynumber(~)
            p = 400;
        end
        
        function b = isAvailable(obj, settings)
            initialpoint = settings.IP;
            b = any(strcmp(initialpoint.getILabel(), obj.pointlabel));
        end

        function msg = check(obj, settings)
            parametersmodel = settings.getSetting('parameters');
            activeParams = parametersmodel.getActive();
            bp = parametersmodel.getBranch();
            if (length(activeParams) ~= obj.nrActive)
                msg = sprintf('You have to select exactly %s free parameters', num2str(obj.nrActive));
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
            paramsetting = settings.getSetting('parameters');
            bp = paramsetting.getBranch();            
            [x0,v0] = obj.initFunction(handle, x0,  param, activeParams, bp);
            v0 = []; %FIXME TODO: make it consistent with call in (old)GUI.
        end
    end
end