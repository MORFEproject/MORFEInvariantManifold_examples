
classdef SimCompSolution < CompSolution
    
    properties
       t
       y
       method
       tspan
       x0
       options
       param
       tE
       yE
       iE        
    end
    properties(Hidden)

       connectData
    end
    
    methods
        function obj = SimCompSolution(settings, compbranch, t, y, tE, yE, iE, method, tspan, x0, options, param)
           obj = obj@CompSolution(settings, compbranch);
           obj.t = t; obj.y = y; obj.method = method; obj.tspan = tspan; obj.x0 = x0; obj.options = options; obj.param = param;
           obj.tE = tE; obj.yE = yE; obj.iE = iE;
           
           obj.connectData = [];
        end
        
        function switches = listSwitches(obj, verbose)
            if nargin < 2
               verbose = false; 
            end
            
           switches = CLCompSwitches();
           if ~isempty(obj.connectData)
                switches.addSwitch(obj.connectData.label, 'Select Connection', @() obj.switchToConnection(), struct('selectmsg', ['Select entire orbit as ' obj.connectData.label], 'interval', diff(obj.tspan), 'coord', obj.connectData.saddle_point  , 'param', obj.param, 'time', obj.t(1),  'toString', @(s) sprintf('---'))); 
           else
                switches.addSwitch('LC', 'Select Cycle', @() obj.switchToLC(), struct('selectmsg', 'Select entire orbit as limit cycle (LC)'  ,'interval', diff(obj.tspan), 'coord', obj.y(1, :), 'param', obj.param, 'time', obj.t(1),  'toString', @(s) sprintf('---'))); 
           end
           
           obj.addSwitches(switches, 1);
           if verbose
                switches = obj.addSwitches(switches, 2:size(obj.y, 1)-1);
           else
                switches = obj.addEventSwitches(switches);
           end
           obj.addSwitches(switches, size(obj.y, 1));
        end
        
        function switches = addEventSwitches(obj, switches)
            for k = 1:length(obj.tE)
                data = struct('npoint', [], 'coord', obj.yE(k, :), 'param', obj.param, 'time', obj.tE(k), 'toString', @(s) sprintf('%g: %s', s.time, mat2str(s.coord)));
                switches.addSwitch(['E',  num2str(obj.iE(k))], ['Event ' num2str(obj.iE(k))], @() obj.switchAtEvent(k), data)
            end
        end
        
        
        function switches = addSwitches(obj, switches, range)
           for index = range
               coord = obj.y(index, :);
               paramvals = obj.param;
               time = obj.t(index);
               
               data = struct('npoint', index, 'coord', coord, 'param', paramvals, 'time', time, 'toString', @(s) sprintf('%g: %s', s.time, mat2str(s.coord)));
               
               if index == 1
                    label = 'P'; msg = 'First Point';
               elseif index == size(obj.y, 1)
                    label = 'P'; msg = 'Last Point';
               else
                    label = 'P'; msg = 'Point';
               end 
               switches.addSwitch(label, msg, @() obj.switchAtPoint(index), data);    
           end
        end
        
        
        function settings = switchAtPoint(obj, index)
            settings = obj.settings.copy();
            settings.setValue('coord', obj.y(index, :));
            settings.setValue('time', obj.t(index));
        end
        
        function settings = switchAtEvent(obj, index)
            settings = obj.settings.copy();
            settings.setValue('coord', obj.yE(index,:));
            settings.setValue('time', obj.tE(index));            
        end
        
        function settings = switchAtEventByTime(obj, time, eventlabel)
            i = str2double(eventlabel(isstrprop(eventlabel,'digit')));
            eventindices = 1:length(obj.tE);
            eventindices = eventindices(obj.iE == i);
            [~, ii] = min(abs(obj.tE(eventindices) - time));
            index = eventindices(ii);
            settings = obj.switchAtEvent(index);
            
        end
        
        
        function settings = switchToLC(obj)
            settings = []; %switch failed
            
            prompt  = {'Enter the tolerance to select a cycle';'Enter the number of test intervals'};
            title   = 'Choose tolerance and ntst';
            lines   = 1;
            def     = {'1e-2', num2str(DefaultValues.STARTERDATA.ntst)};
            answer  = inputdlg(prompt,title,lines,def);           
            if isempty(answer)||isempty(str2double(answer{1}))
                return;
            else
                tolerance = str2double(answer{1});
                ntst=str2double(answer{2});
            end
            x = obj.y';
            t = obj.t;
            
            xstart=round(size(x,2)/3);
            amin=sum(abs(x(:,xstart:end)-x(:,1)*ones(1,size(x(:,xstart:end),2))));
            ep=tolerance;
            amin(find(amin(:)<ep))=inf;
            [pp,qq]=min(amin-(1e-4));
            if (pp==Inf) || pp<0
                ind=[];
            else
                ind = qq(1)+xstart-1;
            end
            if isempty(ind)
                warndlg('No cycle can be found!')
                return;
            end
            
            x=x(:,1:ind);
            t=t(1:ind)';
            tn=(t-t(1))/(t(end)-t(1));
            
            [a,x,tn] = newmeshcycle(x,tn,size(x,2)-1,1,ntst,4);
            x = interp(tn,1,x,a,4);
            sdata = struct();
            sdata.data = struct();
            sdata.data.timemesh = a;
            sdata.data.ntst = ntst;
            sdata.data.ncol = 4;
            sdata.data.parametervalues = obj.param';
            
            x = reshape(x,size(x,2)*size(x,1),1);
            x(end+1) = t(end)-t(1);
            sdata.data.T = x(end);

            % Fill in first parameter, this gets changed later by init_LC_LC when the true parameter selection is made by the user
            x(end+1) = obj.param(1);  %select first parameter as placeholder
      
            sdata.label = 'LC';
            sdata.index = 1;
            sdata.msg = 'LC from Orbit';
            source = struct('x', x, 'v', [], 'restoreGlobals', @() []);
            
            
            settings = obj.settings.copy();
            %add ntst if not exist
            settings.addSetting('ntst', CLSetting('ntst', DefaultValues.STARTERDATA.ntst, InputRestrictions.INT_g0, 2, 4, 1, CLSettingsHelp.getHelp('ntst')));
            %update ntst
            settings.setValue('ntst', ntst);
            settings.setValue('IP', sdata);
            settings.getSetting('IP').setSource(source);
            
            % % % % %
            
            
        end
        function settings = switchToConnection(obj)
            [settings, msg] = obj.connectData.loader(obj);
            if ~isempty(msg)
               fprintf(2, '%s\n', msg); 
            end
            
        end
        
        function save(obj, filename)
            data = struct('t', obj.t, 'y', obj.y, 'tE', obj.tE, 'yE', obj.yE, 'iE', obj.iE, ...
                'method', obj.method, 'tspan', obj.tspan, 'x0', obj.x0, 'options', obj.options, ...
                'param', obj.param, 'connectData', obj.connectData, ...
                'gui', struct('settings', obj.settings,'compbranch', obj.compbranch, 'loader', @(dat) SimCompSolution.loader(dat)));
            
            save(filename, '-struct', 'data')
        end
        function pr = getPreviewRenderer(obj)
            pr = GUIPreviewRendererSimCurve();
        end

        function [b, msg] = showTable(obj, windowname, session)
            if nargin < 3
               session = []; 
            end
            b = 1; msg = '';
            if nargin < 2
                windowname = '';
            end
            GUISimCurveTable.previewPanel(obj, windowname, session)

        end        
        
    end
    methods(Static)
        function obj = loader(dat)
            obj = SimCompSolution(dat.gui.settings, dat.gui.compbranch, dat.t, dat.y, dat.tE, dat.yE, dat.iE, dat.method, dat.tspan, dat.x0, dat.options, dat.param);
            if isfield(dat, 'connectData')
               obj.connectData = dat.connectData; 
            end
        end
    end
end
