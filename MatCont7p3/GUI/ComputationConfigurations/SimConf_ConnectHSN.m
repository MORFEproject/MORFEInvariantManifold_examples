classdef SimConf_ConnectHSN < SimConf_ConnectBase

    
    methods
        
        function obj = SimConf_ConnectHSN(name, type, simconf)
            obj = obj@SimConf_ConnectBase(name, type, simconf, @GUIConnect_getPointHSN);       
            obj.pointlabel = 'LP';
        end
        
        function configureSettings(obj, settings)
            configureSettings@SimConf_ConnectBase(obj, settings);
            settings.addSetting('con_UParam1', CLSetting('UParam1', 0   , InputRestrictions.NUM, 2, 3, 1, ''));
            settings.addSetting('con_eps0'   , CLSetting('eps0'   , 0   , InputRestrictions.POS, 2, 3, 3, ''));
        end
        
    end
    
    
end
