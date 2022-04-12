classdef NumericConfig < CLSettingInterface
    properties
        categories;
        labelstructure;
    end
    events
        visibilityChanged;
        categoriesChanged;
    end
    
    
    methods
        function obj = NumericConfig()
            obj.categories = struct();
            obj.reset();
        end
        
        function declareCategory(obj, label, catname, sortkey, defaultselected)
            if ~isfield(obj.categories, label)
                obj.categories.(label) = struct('sortkey', sortkey, 'valid', 1, 'selected',  defaultselected, 'name', catname);
            else
                obj.categories.(label).valid = 1;
            end
        end
        
        
        
        function reset(obj)
            names = fieldnames(obj.categories);
            for index = 1:length(names)
                obj.categories.(names{index}).valid = 0;
                
            end
            obj.labelstructure = struct();
            
        end
        
        function setLabels(obj,catlabel, labels)
            obj.labelstructure.(catlabel) = labels;
        end
        
        function names = getCategories(obj)
            names = fieldnames(obj.categories);
            f = @(name) obj.categories.(name).valid;
            names = names(find(cellfun(f, names)));
            f = @(name) obj.categories.(name).sortkey;
            [~, ii] = sort(cellfun(f, names));
            names = names(ii);
        end
        
        function b = isCatVisible(obj, catlabel)
            b = obj.categories.(catlabel).selected;
        end
        function  setCatVisible(obj, catlabel, bool)
            obj.categories.(catlabel).selected = bool;
            obj.notify('visibilityChanged');
        end
        
        
        function name = getFullName(obj, catlabel)
            name = obj.categories.(catlabel).name;
        end
        function ls = getLabels(obj, catlabel)
            ls = obj.labelstructure.(catlabel);
        end
        
        function configDone(obj)
            obj.notify('categoriesChanged');
        end
        

        
        % % % CLSetting INTERFACE % % %
        function newobj = copy(obj, ~)
            newobj = NumericConfig();
            newobj.categories = obj.categories;
            newobj.labelstructure = obj.labelstructure;
        end
    end

    
end