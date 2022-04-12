classdef ContConf_NS_ZH < ContConf_NS
    
    methods
        function obj =  ContConf_NS_ZH()
            obj = obj@ContConf_NS('ZH', @init_ZH_NS);
        end
   
         function s = getInitStr(obj)
            s = sprintf('%s, conditional', func2str(obj.initFunction));
        end       
        
        
        function configureSettings(obj, settings)
            configureSettings@ContConf_NS(obj, settings);
            settings.addSetting('eps', CLSetting('amplitude', 1e-6, InputRestrictions.POS, 2, 5, 1, CLSettingsHelp.getHelp('eps') ));

            parammodel = settings.getSetting('parameters');
            parammodel.setEditable(true, true); %values: editable, active: editable
        end
        
    end
    
    methods(Access=protected)
        
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
            
            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
            eps = settings.eps;

            [x0,v0] = obj.initFunction(handle, IP.source.x, param,  sdata, activeParams,  settings.ntst, settings.ncol, eps);
        end
    end
end