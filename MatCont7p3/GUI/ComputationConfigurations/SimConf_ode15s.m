classdef SimConf_ode15s < SimConf
    %continuation configuration
    
    methods
        function obj = SimConf_ode15s(refine, allowEvents, priority)
            obj = obj@SimConf('ode15s', refine, allowEvents, priority);
        end

        function configureSettings(obj, settings)
            configureSettings@SimConf(obj, settings);
            obj.install(settings, 'BDF', false, InputRestrictions.BOOL, [3, 1, 108]);
            obj.install(settings, 'MaxOrder', '5', InputRestrictionCategory({'1', '2', '3', '4', '5'}), [3, 1, 109]);
        end

    end
    methods(Access=protected)
        function [tspan, options] = buildOptions(obj, settings)
            [tspan, options] = buildOptions@SimConf(obj, settings);
            options = odeset(options, 'BDF', CLbool2text(settings.BDF));
            options = odeset(options, 'MaxOrder', str2double(settings.MaxOrder));
            
        end
        
    end
    
end
