classdef ContConf_Hom < ContConf
    properties
        initFunction;
    end
    
 
    methods
        function obj = ContConf_Hom() 
            obj.curvedefinition = @homoclinic;
            obj.label = 'Hom';
            obj.defaultPointType = 'Hom';
            
            %{'NS ';'DRS';'DRU';'NDS';'NDU';'3LS';'3LU';'SH ';'NCH';'BT ';'OFS';'OFU';'IFS';'IFU'};
            obj.testLabels = {'NS',  @(s) ContConf.dimensionCheck(s, 2);...
                              'DRS', @(s) ContConf.dimensionCheck(s, 2);...
                              'DRU', @(s) ContConf.dimensionCheck(s, 2);...
                              'NDS', @(s) ContConf.dimensionCheck(s, 2);...
                              'NDU', @(s) ContConf.dimensionCheck(s, 2);...
                              'test_3LS', @(s) ContConf.dimensionCheck(s, 2);...
                              'test_3LU', @(s) ContConf.dimensionCheck(s, 2);...
                              'SH',  @(s) ContConf.dimensionCheck(s, 2);...
                              'NCH', @(s) ContConf.dimensionCheck(s, 2);...
                               'BT', @(s) ContConf.dimensionCheck(s, 2);...
                              'OFS', @(s) ContConf.dimensionCheck(s, 2);...
                              'OFU', @(s) ContConf.dimensionCheck(s, 2);...
                              'IFS', @(s) ContConf.dimensionCheck(s, 2);...
                              'IFU', @(s) ContConf.dimensionCheck(s, 2)};
                         
            obj.initFunction = @init_Hom_Hom;
            
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf(obj) '_' obj.label];
        end
        
        function s = getInitStr(obj)
            s = func2str(obj.initFunction);
        end
        
        
        function list = getGlobalVars(obj)
            list = {'cds', 'homds'};
            
        end

        function b = isAvailable(~, settings)
            initialpoint = settings.IP;
            plabel = initialpoint.getILabel();
            b =  ~isempty(regexp(lower(plabel), '^hom', 'once'));
        end
        
        function bool = hasPeriod(~)
           bool = 1; 
        end
        
        function configureSettings(obj, settings)
            configureSettings@ContConf(obj, settings);
            obj.setupInitialPoint(settings, false);
            obj.setupTestFunctions(settings, obj.testLabels, false);
            
            % % % Homoclinic parameters
            obj.install(settings, 'ntst', DefaultValues.STARTERDATA.ntst, InputRestrictions.INT_g0, [2, 4, 1]);
            obj.install(settings, 'ncol', DefaultValues.STARTERDATA.ncol, InputRestrictions.INT_g0, [2, 4, 2]);
            
            
            homoclinicparameters = settings.getSetting('homoclinicparameters');
            if isempty(homoclinicparameters)
                homoclinicparameters = CLSettingCustomParameters(settings, {'T', 'eps0', 'eps1'}, true);
                settings.addSetting('homoclinicparameters', homoclinicparameters);
            else
                homoclinicparameters.revive(settings, true);
            end
            
            obj.install(settings, 'eigenvalues', true, InputRestrictions.BOOL, [2, 10, 1]);
            

        end
        
        
        function options = constructContOptions(obj, options, settings)
            options = constructContOptions@ContConf(obj, options, settings);
            options.Eigenvalues = settings.eigenvalues;
        end
        
        function p = getPrioritynumber(~)
            p = 9999;
        end
        
        function msg = check(obj, settings)
            msg = '';
            parametersmodel = settings.getSetting('parameters');
            activeParams = parametersmodel.getActive();
            
            if length(activeParams) ~= 2
                msg = 'Two free system parameters are needed';
            end
        end
        
        
        function oi = getOutputInterpreter(~)
            oi = CLContOutputInterpreterConnection();
        end
        
    end

    methods(Access=protected)
        function [valid, msg] = sanityCheck(obj, handle, x0, param, activeParams, settings)
           msg = obj.check(settings);
           valid = isempty(msg);
        end        
        
        function s = postProcess(~, s)  %Hom specific
            for i = 2:length(s)-1
                s(i).internallabel = ['HOM_' strip(s(i).label)];
            end
        end
        
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
            hparam = settings.getSetting('homoclinicparameters');
            extravec = hparam.parametersfree;
            
            [x0, v0] = obj.initFunction(handle, IP.source.x, IP.source.v, sdata, param, activeParams, settings.ntst, settings.ncol, extravec, settings.custom_T, settings.custom_eps0, settings.custom_eps1);
            
            global homds;
            settings.custom_T.set(homds.T);
            
        end        
    end
    
    
end
