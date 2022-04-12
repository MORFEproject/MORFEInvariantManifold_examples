classdef ContConf_Hom_BT < ContConf_Hom

 
    methods
        function obj = ContConf_Hom_BT() 
            obj = obj@ContConf_Hom();
            obj.initFunction = @init_BT_Hom;
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf_Hom(obj) '_BT'];
        end
        

        function b = isAvailable(~, settings)
            initialpoint = settings.IP;
            plabel = initialpoint.getILabel();
            b =  any(strcmp(plabel, 'BT'));
        end
        
        
        function configureSettings(obj, settings)
            configureSettings@ContConf_Hom(obj, settings);
            hparam = settings.getSetting('homoclinicparameters');
            hparam.revive(settings, true, false);
            settings.addSetting('bt_amplitude', CLSetting('Amplitude', 0.0001, InputRestrictions.POS, 2, 5, 1, CLSettingsHelp.getHelp('bt_amplitude') ));
            settings.addSetting('bt_ttolerance', CLSetting('TTolerance', 0.00001, InputRestrictions.POS, 2, 5, 2, CLSettingsHelp.getHelp('bt_ttolerance') ));
        end
        
        function p = getPrioritynumber(~)
            p = 1000000;
        end
        
    end

    methods(Access=protected)     
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
            hparam = settings.getSetting('homoclinicparameters');
            extravec = hparam.parametersfree;
            
            [x0, v0] = obj.initFunction(handle, x0, sdata, param, activeParams, settings.ntst, settings.ncol, settings.bt_ttolerance, settings.bt_amplitude, extravec);

            global homds
            settings.custom_T.set(homds.T);
            settings.custom_eps0.set(homds.eps0);
            settings.custom_eps1.set(homds.eps1);
            
            
        end
    end
    
    
end
