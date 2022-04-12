classdef CLSettingHTDS < CLSettingInterface
    
    properties
       htds = struct(); 
       visible = 0;
       
       setting_changed;
    end
    
    methods
        % % % CLSetting INTERFACE % % %
        function value = getValue(obj)
            value = obj.htds;
        end
        function [valid, msg] = setValue(obj, val)
            valid = 1; 
            msg = '';
            obj.htds = val;
            
            obj.initHTDS();
        end
        function id = getGroupID(~); id = 2; end  %starter
        function id = getSubGroupID(~); id = 3; end
        function id = getItemID(~); id = 999; end
        function s = toString(~); s = 'Homotopy_data'; end
        function h = getHelpStr(~); h = 'contains data used in the homotopy algoritmes'; end
        function b = isVisible(obj); b = obj.visible; end
        function setVisible(obj, v); obj.visible = v; end
        
        function newobj = copy(obj, ~)
            newobj = CLSettingHTDS();
            newobj.htds = obj.htds;
            
        end

        function setItem(obj, name, value)
            value = logical(value);
            if strcmp(name, 'T')
                obj.htds.gds.extravec(1) = value;
            elseif strcmp(name, 'eps1')
                obj.htds.gds.extravec(3) = value;
            else
                index = str2double(name(~isletter(name)));
                if startsWith(name, 'UParam')
                    obj.htds.gds.UParamsFree(index) = value;
                elseif startsWith(name, 'SParam')
                    obj.htds.gds.SParamsFree(index) = value;
                    
                else
                    assert(false, ['unknown keyword: ', name]);
                end
            end
            
            
            obj.setting_changed();
        end
        
        function b = getItem(obj, name)
            if strcmp(name, 'T')
                b = obj.htds.gds.extravec(1);
            elseif strcmp(name, 'eps1')
                b = obj.htds.gds.extravec(3);
            else
                index = str2double(name(~isletter(name)));
                if startsWith(name, 'UParam')
                    b = obj.htds.gds.UParamsFree(index);
                elseif startsWith(name, 'SParam')
                    b = obj.htds.gds.SParamsFree(index);
                    
                else
                    assert(false, ['unknown keyword: ', name]);
                end
            end
            
            
        end
        function setUParam(obj, index, value)
           obj.htds.gds.UParams{index, 2} = value; 
        end
         function setSParam(obj, index, value)
           obj.htds.gds.SParams{index, 2} = value; 
         end       
         function setT(obj, value)
                obj.htds.T = value;
         end
          function seteps1(obj, value)
                obj.htds.eps1 = value;
         end        


        % render UParam/Sparam/T/eps selections
        function box = renderGUI(obj, session, settings, label, panelhandle, options, suggestions)
            
            obj.setting_changed = @() settings.refresh();
            
            
            
            items = [obj.htds.gds.UParams; obj.htds.gds.SParams; {'T',  obj.htds.T}; {'eps1', obj.htds.eps1}];
            
            getters = {};
            setters = {};
            for k = 1:size(obj.htds.gds.UParams, 1)
                getters{end+1} = @() obj.htds.gds.UParams{k, 2};
                setters{end+1} = @(val) obj.setUParam(k, val);
            end
             for k = 1:size(obj.htds.gds.SParams, 1)
                getters{end+1} = @() obj.htds.gds.SParams{k, 2};
                setters{end+1} = @(val) obj.setSParam(k, val);
            end           
            getters{end+1} = @() obj.htds.T;
            setters{end+1} = @(val) obj.setT(val);
            getters{end+1} = @() obj.htds.eps1;
            setters{end+1} = @(val) obj.seteps1(val);
            grid = cell(size(items, 1), 2);
            
            
            for k = 1:size(items, 1)
                grid{k, 1} = uicontrol(panelhandle, 'Style', 'radiobutton' , 'String', items{k, 1}, 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, 'Value', obj.getItem(items{k, 1}) , 'Callback', @(src, ev) obj.setItem(items{k, 1}, src.Value)  );
                grid{k, 2} = GUIEditBox(panelhandle, getters{k}, setters{k}, InputRestrictions.NUM, 'enable', 'inactive');
                
            end

            
            
            
            box1 = LayoutBox(grid(1:end-2, :));
            box1.widths(1) = max(suggestions.labelsize, box1.widths(1));
            
            box2 = LayoutBox(grid(end-1:end, :));
            box2.widths(1) = max(suggestions.labelsize, box2.widths(1));
      
            sectionsettings = DefaultValues.SECTIONNAMESETTING;
            curvelabel = session.getCompConf().getCurveLabel();
      
            if contains(curvelabel, 'Het')
               header = 'Heteroclinic parameters';
            else
                header = 'Homoclinic parameters';
            end
            
            
            grid = {GUIEmptySpace(round(DefaultValues.LETTERDIMENSION(1)), 10); ...
                uicontrol(panelhandle, 'Style' , 'Text', 'String', 'Connection parameters', sectionsettings{:});...
                box1; ...
                %GUIEmptySpace(round(DefaultValues.LETTERDIMENSION(1)/2), 10); ...
                uicontrol(panelhandle, 'Style' , 'Text', 'String', header, sectionsettings{:}); ...
                box2; ...
                uicontrol(panelhandle, 'Style' , 'Text', 'String', 'eps1 tolerance', sectionsettings{:})};
                 
            
            
            
            box = LayoutBox(grid);
        end
        
        
        function initHTDS(obj)
            assert(isfield(obj.htds.gds, 'SParams'), 'unexpected: SParams is missing');
            assert(isfield(obj.htds.gds, 'UParams'), 'unexpected: UParams is missing');
            assert(isfield(obj.htds.gds, 'extravec'), 'unexpected: extravec is missing');
            assert(isfield(obj.htds, 'T'), 'unexpected: T is missing');
            assert(isfield(obj.htds, 'eps1'), 'unexpected: eps1 is missing');
            
            if ~isfield(obj.htds.gds, 'SParamsFree') || length(obj.htds.gds.SParamsFree) ~= size(obj.htds.gds.SParams, 1)   
               obj.htds.gds.SParamsFree =  repmat([false], 1, size(obj.htds.gds.SParams, 1) );
            end
            if ~isfield(obj.htds.gds, 'UParamsFree') || length(obj.htds.gds.UParamsFree) ~= size(obj.htds.gds.UParams, 1)   
               obj.htds.gds.UParamsFree =  repmat([false], 1, size(obj.htds.gds.UParams, 1) );
            end           
            if ~isfield(obj.htds.gds, 'period')
               obj.htds.gds.period = 1; 
            end
            
        end
        
        function cleanupSParam(obj, settings)

            dim_nneg = size(obj.htds.gds.SParams, 1);
            SParamTestTolerance = settings.getVal('SParamTestTolerance' , 1e-05);
            
            %{
            %NN/WG:  eliminated, use unknown. Does not know why this is
            %here
            
            for i = 1:dim_nneg
                prod = 1;
                stablepars = [];
                index = [];
                for i = 1:dim_nneg
                    if abs(obj.htds.gds.SParams{i,2})~=0
                        prod = prod*abs(obj.htds.gds.SParams{i,2});
                        stablepars = [stablepars;abs(obj.htds.gds.SParams{i,2})];
                        index = [index;i];
                    end
                end
                if prod < SParamTestTolerance
                    [~,pos] = min(stablepars);
                    obj.htds.gds.SParams{index(pos),2} = 0;
                end
            end
            %}
            
            for i = 1:dim_nneg
                if abs(obj.htds.gds.SParams{i,2}) < SParamTestTolerance
                    obj.htds.gds.SParams{i,2} = 0;
                end
            end
            
            obj.htds.gds.extravec = obj.htds.gds.extravec * 0;
            obj.htds.gds.SParamsFree = obj.htds.gds.SParamsFree * 0;
            obj.htds.gds.UParamsFree = obj.htds.gds.UParamsFree * 0;
        end
        
        
    end
end
