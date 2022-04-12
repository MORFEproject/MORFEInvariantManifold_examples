classdef ContConf_LCCodim2 < ContConf
    properties
        initFunction;
    end
    
 
    methods

        function s = getLabel(obj)
            s = [getLabel@ContConf(obj) '_' obj.label];
        end
        
        function s = getInitStr(obj)
            s = func2str(obj.initFunction);
        end
        
        function configureSettings(obj, settings)
            configureSettings@ContConf(obj, settings);
            
            obj.setupInitialPoint(settings, false); %false: don't show coords by default.
            obj.setupTestFunctions(settings, obj.testLabels);
            obj.install(settings, 'multipliers', true, InputRestrictions.BOOL, [2, 9, 1]);

            %Add period, but period not editable
            settings.addSetting('Period', CLSettingPeriod(false));
            o = settings.getSetting('Period'); o.setSelectable(false);
            
            %add discretization
            obj.install(settings, 'ntst', DefaultValues.STARTERDATA.ntst, InputRestrictions.INT_g0, [2, 4, 1]);
            obj.install(settings, 'ncol', DefaultValues.STARTERDATA.ncol, InputRestrictions.INT_g0, [2, 4, 2]);

        end
        
        function options = constructContOptions(obj, options, settings)
            %settings to cont-struct:
            options = constructContOptions@ContConf(obj, options, settings); %superclass
            options.Multipliers = settings.multipliers; %add multipliers
            
        end
        
        function p = getPrioritynumber(~)
            p = 1000;
        end
        
        function bool = hasPeriod(~)
           bool = 1; 
        end
        function msg = check(obj, settings)
            msg = '';
            parametersmodel = settings.getSetting('parameters');
            activeParams = parametersmodel.getActive();
            if length(activeParams) ~= 2
                msg = 'You have to select two free parameters';
            end
        end
    end
    methods(Access=protected)
        function [valid, msg] = sanityCheck(obj, handle, x0, param, activeParams, settings)
            valid = true;
            msg = '';
            if length(activeParams) ~= 2
                valid = false;
                msg = 'You have to select two free parameters';
            end
        end        
        
    end
    
    
end
