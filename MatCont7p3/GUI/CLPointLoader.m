classdef CLPointLoader < handle
    properties
        currentlabel;
        labelmap;
        session;
        
        clicktimestamp = 0;
    end
    
    
    methods
        function obj = CLPointLoader(session)
            if nargin < 1
                session = [];
            end
            obj.session = session;
            obj.currentlabel = '';
            obj.clear();
        end
        
        function clear(obj)
            obj.labelmap = containers.Map();
            obj.clicktimestamp = clock();
        end
        
        function newLabel(obj)
            obj.currentlabel = sprintf('lbl%i', tic());
        end

        function l = getCurrentLabel(obj)
            l = obj.currentlabel;
        end
        
        function l = getLabel(obj, path)
            keys = obj.labelmap.keys();
            for i = 1:length(keys)
                if strcmp(path, obj.labelmap(keys{i}))
                    l = keys{i};
                    return;
                end
            end
            
            l = sprintf('lbl%i', tic());
            obj.labelmap(l) = path;
            
        end
        
        function assignPath(obj,  filepath)
            obj.labelmap(obj.getCurrentLabel()) = filepath;
        end        
        
        function renaming(obj, oldfilepath, newfilepath)
            keys = obj.labelmap.keys();
            for i = 1:length(keys)
                if strcmp(oldfilepath, obj.labelmap(keys{i}))
                    obj.labelmap(keys{i}) = newfilepath;
                    return;
                end
            end
        end
        
        function disp(obj)
            keys = obj.labelmap.keys();
            fprintf('current: %s\n\n', obj.currentlabel);
            for i = 1:length(keys)
                fprintf('%s   %s\n', pad(keys{i}, 20), obj.labelmap(keys{i}));
            end
        end
        
        function loadPoint(obj, label, sdata)
            if isempty(obj.session); return; end
                        
            
            if ~isvalid(obj.session)
                fprintf(2, 'Connection with current session has been lost, unable to load in point\n');
                return;
            end
            
            if obj.labelmap.isKey(label)
                path = obj.labelmap(label);
                %assumed that solution is continuation curve.
                
                if exist(path, 'file') == 2 %file exists and is matlab
                    
                    timestamp = clock();
                    timebetweenclicks = etime(timestamp, obj.clicktimestamp);
                    obj.clicktimestamp = timestamp;
                    
                    pointlabel = strip(sdata.label);
                    [~, filename, ~] = fileparts(path);
                    if timebetweenclicks < 0.6
                        
                        
                        if strcmp('Yes', questdlg(sprintf('Load in %s (%s) from "%s"?', sdata.msg, pointlabel, filename), sprintf('Select %s as a new initial point', pointlabel)))
                            
                            contcurve = CompSolution.load(path);
                            newsettings = contcurve.switchAtPoint(sdata.index);
                            obj.session.changeState(obj.session.getSystem(), obj.session.getDiagram(), contcurve, newsettings);
                        end
                    else
                        content = GUIPreviewRendererContCurve.getBifurcationData(sdata, @toString);
                        content = [{'npoint', num2str(sdata.index)}; content];
                        
                        labellen =  max(cellfun(@(x) length(x), content(:, 1)));
                        
                        fprintf('*** %s: %s (%s)\n',  filename, sdata.msg, pointlabel);
                        for k = 1:size(content, 1)
                            fprintf('         %s %s\n', pad([content{k, 1} ':'], labellen + 2), content{k, 2});
                        end
                    end
                    return;
                end
            end
            fprintf(2, 'Curve can no longer be found\n')
            
        end
        function loadPointSim(obj, label, data)
            if isempty(obj.session); return; end
                        
            
            if ~isvalid(obj.session)
                fprintf(2, 'Connection with current session has been lost, unable to load in point\n');
                return;
            end
            
            if obj.labelmap.isKey(label)
                path = obj.labelmap(label);
                %assumed that solution is simulation curve.
                
                if exist(path, 'file') == 2 %file exists and is matlab
                    
                    timestamp = clock();
                    timebetweenclicks = etime(timestamp, obj.clicktimestamp);
                    obj.clicktimestamp = timestamp;
                    
                    t = data{1}; y = data{2}; eventlabel = data{3};
                    
           
                    [~, filename, ~] = fileparts(path);
                    if timebetweenclicks < 0.6
                        
                        
                        if strcmp('Yes', questdlg(sprintf('Load in %s (at t=%6.6g) from "%s"?', eventlabel, t, filename), sprintf('Select %s as a new initial point', eventlabel)))
                            
                            simcurve = CompSolution.load(path);
                            newsettings = simcurve.switchAtEventByTime(t, eventlabel);
                            obj.session.changeState(obj.session.getSystem(), obj.session.getDiagram(), simcurve, newsettings);
                            
                        end
                    else
                        fprintf('%s: t=%5.5g, y=[%s]\n', eventlabel, t, sprintf('%5.5g ', y));
                    end
                    return;
                end
            end
            fprintf(2, 'Curve can no longer be found\n')
            
        end            
            
     
        
    end
end

function s = toString(x)
    if ischar(x)
        s = x;
    else
       s =  sprintf('%.7e ',x);
    end
end
