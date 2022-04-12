classdef SimConf_ConnectHet < SimConf_ConnectBase

    
    methods
        
        function obj = SimConf_ConnectHet(name, type, simconf)
            obj = obj@SimConf_ConnectBase(name, type, simconf, @GUIConnect_getPointHet);            
        end
        
        function configureSettings(obj, settings)
            configureSettings@SimConf_ConnectBase(obj, settings);
            settings.addSetting('con_UParam1', CLSetting('UParam1', 0   , InputRestrictions.NUM, 2, 3, 1, ''));
            settings.addSetting('con_UParam2', CLSetting('UParam2', 0   , InputRestrictions.NUM, 2, 3, 2, ''));
            settings.addSetting('con_eps0'   , CLSetting('eps0'   , 0   , InputRestrictions.POS, 2, 3, 3, ''));
            
            
            system = settings.system;
            coordinates = CLSettingCustomCoordinates(settings, system.getCoordinates(), 'coord_target', 'co2_', 131);
            settings.addSetting('coord_target', coordinates);          

        end
        
    end
    
    
end
