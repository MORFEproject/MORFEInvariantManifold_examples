classdef ContConf_Hom_Hom < ContConf_Hom

 
    methods

        
        function configureSettings(obj, settings)
            configureSettings@ContConf_Hom(obj, settings);
            initial = settings.getSetting('IP');
            newT = initial.currentpointdata.data.T;
            settings.custom_T.set(newT);
            settings.custom_eps0.set(initial.currentpointdata.data.meta.eps0);
            settings.custom_eps1.set(initial.currentpointdata.data.meta.eps1);
            
        end
        
    end
    
    
end
