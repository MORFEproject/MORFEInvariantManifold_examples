classdef ContConf_LPC_GH < ContConf_LPC
    
    methods
        function obj =  ContConf_LPC_GH()
            obj = obj@ContConf_LPC('GH', @init_GH_LPC);
            
            
        end
        
        function configureSettings(obj, settings)
            configureSettings@ContConf_LPC(obj, settings);
            
            settings.addSetting('eps', CLSetting('amplitude', 1e-6, InputRestrictions.POS, 2, 5, 1, CLSettingsHelp.getHelp('eps') ));
            %fprintf(2, 'TODO: re-enable parameter values when selecting \n');
              
            parammodel = settings.getSetting('parameters');
            parammodel.setEditable(true, true); %values: uneditable, active: editable
                 
            
            
        end
        
        
    end
    
        methods(Access=protected)
        
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
            paramsetting = settings.getSetting('parameters');
            bp = paramsetting.getBranch();            
            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
            
            eps = settings.eps;
            
            [x0,v0] = obj.initFunction(handle, IP.source.x, param,  sdata, activeParams,  settings.ntst, settings.ncol, eps, bp);
        end
    end
    
    
end
