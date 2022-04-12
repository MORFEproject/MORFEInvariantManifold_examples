function solution = LegacySolutionLoader(filepath)
            [syspath, ~, ~] = fileparts(filepath);
            [syspath, ~, ~] = fileparts(syspath);
            system = CLSystem([syspath '.mat']);
            
        
            data = load(filepath);
            pointtype =  strip(data.point);
            curvetype = strip(data.ctype);
            
            
            COMPCONFIGS = CompConfigurations();
            k = 1;
            while k <= length(COMPCONFIGS.conflist) && ~strcmp(curvetype, COMPCONFIGS.conflist{k}.getSolutionLabel())
                k = k + 1;
            end
            if k > length(COMPCONFIGS.conflist)
               fprintf(2, 'unknown curvetype: %s\n', curvetype); 
               solution = [];
               return;
            end
            
            compbranch = COMPCONFIGS.conflist{k};
          
            
            settings = CLSettings();
            settings.installSetting('system', system);
            compbranch.configureSettings(settings);
            settings.setValue('IP', pointtype);
            
            
            P0 = [];
            ap = [];
            names = fieldnames(data);
            for k = 1:length(names)
                name = names{k};
                if regexp(name, 'ds$')
                    eval(sprintf('global %s; %s = data.%s;', name, name, name));
                    if isfield(data.(name), 'P0')
                        P0 = data.(name).P0;
                        if isfield(data.(name),'ActiveParams')
                           ap =  data.(name).ActiveParams;
                        end
                    end
                end
            end
            if ~isempty(P0)
                settings.setValue('parameters', P0);
            end
            if ~isempty(ap)
               a = settings.getSetting('parameters');
               a.parametersfree(ap) = 1;
            end
            

            
            
            if strcmp(curvetype, 'O')
                settings.setValue('coord', data.x(:,1));
                settings.setValue('parameters', data.param);
                solution = SimCompSolution(settings, compbranch, data.t, data.x',[], [], str2func(compbranch.getMethodName()), [data.t(1), data.t(end)], data.x(:,1), data.option, data.param);

            else
                solution = ContCurve(settings, compbranch, data.x, data.v, data.s, data.h, data.f);
                oi = compbranch.getOutputInterpreter();
                [coord, param] = oi.getPoint(solution, 1);
                solution.settings.setValue('coord', coord);
                solution.settings.setValue('parameters', param);
            end

    
end
