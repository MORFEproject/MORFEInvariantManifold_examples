classdef ContConf_LC < ContConf
    properties
        initFunction;
    end
    
 
    methods
        function obj = ContConf_LC() 
            obj.curvedefinition = @limitcycle;
            obj.label = 'LC';
            obj.defaultPointType = 'LC';
            
            obj.testLabels = {'BPC', @(s) ContConf.dimensionCheck(s, 1); 'PD', @(s) ContConf.dimensionCheck(s, 1); 'LPC', @(s) ContConf.dimensionCheck(s, 1); 'NS', @(s) ContConf.dimensionCheck(s, 2)};
                
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf(obj) '_' obj.label];
        end
        
        function s = getInitStr(obj)
            s = func2str(obj.initFunction);
        end
        
        
        function list = getGlobalVars(obj)
            list = {'cds', 'lds'};
            
        end

        function b = isAvailable(~, settings)
            b = 0;
        end
        function bool = hasPeriod(~)
           bool = 0; 
        end
        function configureSettings(obj, settings)
            configureSettings@ContConf(obj, settings);
            
            
            obj.setupTestFunctions(settings, obj.testLabels);
            obj.install(settings, 'multipliers', true, InputRestrictions.BOOL, [2, 9, 1]);
            
            settings.addSetting('Period', CLSettingPeriod(true));
            o = settings.getSetting('Period'); o.setSelectable(true);
            obj.install(settings, 'ntst', DefaultValues.STARTERDATA.ntst, InputRestrictions.INT_g0, [2, 4, 1]);
            obj.install(settings, 'ncol', DefaultValues.STARTERDATA.ncol, InputRestrictions.INT_g0, [2, 4, 2]);
            system = settings.system;
            
            dim = length(system.getCoordinates());
            
            %FIXME TODO: prcInput is system afhankelijk, dus moet
            %verwijderd worden bij systemdimension change.
            settings.addSetting('prcInput', CLSetting('Input', 1, InputRestrictions.vectorOrScalar(dim), 2, 100, 2, ''));
            settings.addSetting('PRCenabled', CLSettingPRCWindow('PRC'));
            settings.addSetting('dPRCenabled', CLSettingPRCWindow('dPRC'));
            
            s = settings.getSetting('PRCenabled'); s.syncValue();
            s = settings.getSetting('dPRCenabled'); s.syncValue();
            
        end
        function options = constructContOptions(obj, options, settings)
            options = constructContOptions@ContConf(obj, options, settings);
            options.Multipliers = settings.multipliers;
                        
            options.PRC = settings.PRCenabled;
            options.dPRC = settings.dPRCenabled;
            options.Input = settings.prcInput;
        end
        
        function p = getPrioritynumber(~)
            p = 1000;
        end
        
        function msg = check(obj, settings)
            msg = '';
            parametersmodel = settings.getSetting('parameters');
            activeParams = parametersmodel.getActive();
            if length(activeParams) + settings.Period ~= 2
                msg = 'You have to select either one free parameter and the period or two free parameters';
            end
        end
        
        
    end

    methods(Access=protected)
        function [valid, msg] = sanityCheck(obj, handle, x0, param, activeParams, settings)
            valid = true;
            msg = '';
            if length(activeParams) + settings.Period ~= 2
                valid = false;
                msg = 'You have to select either one free parameter and the period or two free parameters';
            end
        end
        
        function s = postProcess(~, s)  %LC specific
            for i = 1:length(s)
                if strcmp(strip(s(i).label), 'NS') && contains(lower(s(i).msg), 'neutral saddle')
                    s(i).internallabel = 'NC';
                end
            end
        end
    end
    
    
end
