classdef GUILineOptions < CLSettingInterface
    properties(Constant, Hidden)
        DEFAULT_LINE_SETTING = {'Color', 'blue', 'linestyle', '-'};
        DEFAULT_MARKER_SETTING = {'Color', 'red', 'marker', '*'};
        DEFAULT_SPECIAL_SETTING = {'Color', 'red', 'linestyle', '-', 'linewidth' , 3};
        DEFAULT_LABEL_SETTING = {'Color', 'black', 'FontSize', 10};
        
    end
    
    properties
        linesettings;
        defaultlinesettings;
        namemap;

    end
    
    methods
        function obj = GUILineOptions(curvelabels)
            obj.linesettings = struct();
            obj.namemap = containers.Map();
            if nargin > 0
                for k = 1:length(curvelabels)
                    label = curvelabels{k};
                    obj.linesettings.(['CURVE_' label]) = obj.DEFAULT_LINE_SETTING;
                    obj.namemap(['CURVE_' label]) = label;
                    
                    title =  GUICurveModifications.hasMod(label);
                    if ~isempty(title)
                        obj.linesettings.(['CURVE_' label '_MOD']) = {};
                        obj.namemap(['CURVE_' label '_MOD']) = [label ' ' title];
                    end
                    
                    
                end
            end
            obj.linesettings.OTHER_SINGPOINT =  obj.DEFAULT_MARKER_SETTING;
            obj.namemap('OTHER_SINGPOINT') = 'Singularity (Point)';
            
            obj.linesettings.OTHER_SINGCURVE =  obj.DEFAULT_SPECIAL_SETTING;
            obj.namemap('OTHER_SINGCURVE') = 'Singularity (Curve)';
            
            obj.linesettings.OTHER_LABEL=  obj.DEFAULT_LABEL_SETTING;
            obj.namemap('OTHER_LABEL') = 'Singularity (Label)';
            
            obj.defaultlinesettings = obj.linesettings;
            
        end
        
        function o = getLineSetting(obj, curvelabel)
            label = ['CURVE_' curvelabel];
            if isfield(obj.linesettings, label)
                o = obj.linesettings.(label);
            else
                o = obj.DEFAULT_LINE_SETTING;
            end
            
        end
        
        function o = getLineModSetting(obj, curvelabel)
            label = ['CURVE_' curvelabel '_MOD'];
             if isfield(obj.linesettings, label)
                o = obj.linesettings.(label);
            else
                o = {};
            end           
        end
        function o = getSingPointSetting(obj)
           o = obj.linesettings.OTHER_SINGPOINT;
        end
        function o = getSingCurveSetting(obj)
           o = obj.linesettings.OTHER_SINGCURVE;
        end
        function o = getSingLabelSetting(obj)
           o =  obj.linesettings.OTHER_LABEL;
        end
        
        function [o, name] = getSetting(obj, label)
            o = obj.linesettings.(label);
            name = obj.namemap(label);
            
        end
        function setDefault(obj, label)
            obj.linesettings.(label) = obj.defaultlinesettings.(label);
        end
        
        function setSetting(obj, label, options)
            obj.linesettings.(label) = options;
        end
        function l = getLabels(obj)
            l = fieldnames(obj.linesettings);
            
        end
        
        function l = length(obj)
            l = length(fieldnames(obj.linesettings));
        end
        
        %clsettings interface
        function newobj = copy(obj, ~)
            newobj = GUILineOptions();
            newobj.linesettings = obj.linesettings;
            newobj.defaultlinesettings = obj.defaultlinesettings;
            newobj.namemap = obj.namemap;
        end
    end
    
    
    
    
    
    
    
end