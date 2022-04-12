classdef ContConf_NS_HH < ContConf_NS
    
    methods
        function obj =  ContConf_NS_HH()
            obj = obj@ContConf_NS('HH', @init_HH_NS);
        end
        
        function configureSettings(obj, settings)
            configureSettings@ContConf_NS(obj, settings); %configure NS

            settings.addSetting('whichNS', CLSetting('whichNS', '1', InputRestrictionCategory({'1', '2'}), 2, 5, 2, CLSettingsHelp.getHelp('whichNS') ));
            %amplitude is called 'eps' internally
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
            
            if strcmp(settings.whichNS, '1')
                init = @init_HH_NS1;
            else
                init = @init_HH_NS2;
            end
            
            
            [x0,v0] = init(handle, IP.source.x, param,  sdata, activeParams,  settings.ntst, settings.ncol, eps);
        end
    end
end
