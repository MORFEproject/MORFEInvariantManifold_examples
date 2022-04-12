classdef ContConf_LPC < ContConf_LCCodim2
    properties
        pointlabel
        
    end
    
    methods
        function obj =  ContConf_LPC(pointlabel, initFunction)
            obj.curvedefinition = @limitpointcycle;
            obj.label = 'LPC';
            obj.defaultPointType = 'LPC';
            obj.pointlabel = pointlabel;
            obj.initFunction = initFunction;
            obj.testLabels = {'R1', @(s) ContConf.dimensionCheck(s, 2); 'CPC', @(s) ContConf.dimensionCheck(s, 1); 'LPNS', @(s) ContConf.dimensionCheck(s, 3); 'LPPD', @(s) ContConf.dimensionCheck(s, 2)};
            
        end
        function list = getGlobalVars(obj)
            list = {'cds', 'lds'};
            
        end
        
        
        function configureSettings(obj, settings)
            configureSettings@ContConf_LCCodim2(obj, settings);
             paramsetting = settings.getSetting('parameters');
            paramsetting.revive(settings, true, true); %active:true, branch:true    
            
            parammodel = settings.getSetting('parameters');
            parammodel.setEditable(false, true); %values: uneditable, active: editable
      
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf_LCCodim2(obj) '_' obj.pointlabel];
        end
        
        function p = getPrioritynumber(~)
            p = 1;
        end
        
        function b = isAvailable(obj, settings)
            initialpoint = settings.IP;
            b = any(strcmp(initialpoint.getILabel(), obj.pointlabel));
        end
        
    end
    
    methods(Access=protected)
        
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
            paramsetting = settings.getSetting('parameters');
            bp = paramsetting.getBranch();            
            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
            
            
            [x0,v0] = obj.initFunction(handle, IP.source.x,  sdata, activeParams,  settings.ntst, settings.ncol, bp);
        end
    end
end