classdef ContConf_PD < ContConf_LCCodim2
    properties
        pointlabel
        
    end
    
    methods
        function obj =  ContConf_PD(pointlabel, initFunction)
            obj.curvedefinition = @perioddoubling;
            obj.label = 'PD';
            obj.defaultPointType = 'PD';
            obj.pointlabel = pointlabel;
            obj.initFunction = initFunction;
            
            obj.testLabels = {'R2', @(s) ContConf.dimensionCheck(s, 2); 'LPPD', @(s) ContConf.dimensionCheck(s, 2); 'GPD', @(s) ContConf.dimensionCheck(s, 1); 'PDNS', @(s) ContConf.dimensionCheck(s, 3)};
        end
        function list = getGlobalVars(obj)
            list = {'cds', 'lds'};
            
        end
        function s = getLabel(obj)
            s = [getLabel@ContConf_LCCodim2(obj) '_' obj.pointlabel];
        end
        
        
        function p = getPrioritynumber(~)
            p = 2;
        end
        
        function b = isAvailable(obj, settings)
            initialpoint = settings.IP;
            b = any(strcmp(initialpoint.getILabel(), obj.pointlabel));
        end
        function configureSettings(obj, settings)
            configureSettings@ContConf_LCCodim2(obj, settings);
            parammodel = settings.getSetting('parameters');
            parammodel.setEditable(false, true); %values: uneditable, active: editable
            
        end
    end
    
    methods(Access=protected)
        
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
        
            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
  
            [x0,v0] = obj.initFunction(handle, IP.source.x,  sdata, activeParams,  settings.ntst, settings.ncol);
        end
    end
end