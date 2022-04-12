classdef PlotSelections < handle
    
    properties
        categories;
        solutionlabel;
    end
    
    
    methods
        function obj = PlotSelections()
            obj.categories = struct();
        end
        
        
        function declareCategory(obj, catlabel, catname, sortkey, singleton)
           if nargin < 5
              singleton = 0; 
           end
            
           obj.categories.(catlabel) = struct('key', sortkey, 'name', catname, 'options', struct(), 'singleton', singleton); 
            
        end
        
        function declareSubCategory(obj, catlabel, subcatlabel, func)
            obj.categories.(catlabel).options.(subcatlabel) = struct('name',  subcatlabel, 'options', struct());
            if nargin == 4
               obj.declareItem(catlabel, subcatlabel, 'Default', subcatlabel, func); 
            end
        end
        function declareSubCategoryAltName(obj, catlabel, subcatlabel, displayname, func)
            obj.categories.(catlabel).options.(subcatlabel) = struct('name',  displayname, 'options', struct());
            if nargin == 5
               obj.declareItem(catlabel, subcatlabel, 'Default', subcatlabel, func); 
            end
        end
        
        function declareItem(obj, catlabel, subcatlabel, itemlabel, subcatname, func, pointplot)
            if nargin < 7
                pointplot = false;
            end
            obj.categories.(catlabel).options.(subcatlabel).options.(itemlabel) = struct('name', itemlabel, 'subcatname', obj.categories.(catlabel).options.(subcatlabel).name , 'func', func, 'pointplot', pointplot); 
        end
        
        
        function disp(obj)
            catnames = fieldnames(obj.categories);      
            [~, ii] = sort(cellfun( @(x) obj.categories.(x).key    , catnames));
            catnames = catnames(ii);
            
            
            for i = 1:length(catnames)
                catname = catnames{i};
                
                subcatnames = fieldnames(obj.categories.(catname).options);
                
                for j = 1:length(subcatnames)
                    subcatname = subcatnames{j};
                    items = fieldnames(obj.categories.(catname).options.(subcatname).options);
                    for k = 1:length(items)
                        item = items{k};
                        n = obj.categories.(catname).options.(subcatname).options.(item).subcatname;
                        
                        fprintf('%s_%s_%s -> %s\n', catname, subcatname, item, n);
                        
                    end
                end
            end
            fprintf('\n');
        end
        
        
        function [catnames] = getCategories(obj)
            catnames = fieldnames(obj.categories);      
            [~, ii] = sort(cellfun( @(x) obj.categories.(x).key    , catnames));
            catnames = catnames(ii);
            
            for k = 1:length(catnames)
               catnames{k, 2} = obj.categories.(catnames{k}).name;
            end
        end
        
        function [subcatnames, singleton] = getSubCategories(obj, catname)
            subcatnames = fieldnames(obj.categories.(catname).options);
            singleton = obj.categories.(catname).singleton;
        end
        
        function items = getItems(obj, catname, subcatname)
            items = fieldnames(obj.categories.(catname).options.(subcatname).options);
        end
        
        function [itemlist, singleton] = getSubcatItems(obj, catname, itemname)
            itemlist = cell(0, 2);
            subcatnames = obj.getSubCategories(catname);
            singleton = obj.categories.(catname).singleton;
            for k = 1:length(subcatnames)
                subcatname = subcatnames{k};
                items = fieldnames(obj.categories.(catname).options.(subcatname).options);
                for j = 1:length(items)
                    if strcmp(items{j}, itemname)
                        itemlist(end+1, :) = {subcatname, obj.categories.(catname).options.(subcatname).options.(items{j}).subcatname};
                    end
                end
            end
            
        end
        function label = getSelectionLabel(obj, catname, subcatname, item)
           label =  obj.categories.(catname).options.(subcatname).options.(item).subcatname;
        end
        
        function [fhandle, pointplot] = getSelectFunction(obj, catname, subcatname, item)
           fhandle =  obj.categories.(catname).options.(subcatname).options.(item).func;
           pointplot =  obj.categories.(catname).options.(subcatname).options.(item).pointplot;
        end
        
        function b = isValid(obj, catname, subcatname, itemname)
           b = isfield(obj.categories, catname) && isfield(obj.categories.(catname).options, subcatname) && isfield( obj.categories.(catname).options.(subcatname).options, itemname);
        end
        
        function suggestions = getDefaultSelection(obj, session)
            param = session.settings.getSetting('parameters');
            suggestions = cell(0, 3);
            if isfield(obj.categories, 'parameters') && ~isempty(param)
                actives = param.parameters(param.hasActive());
                for k = 1:length(actives)
                   active = actives{k};
                   if isfield(obj.categories.parameters.options, active)
                      suggestions(end+1, :) = {'parameters', active, 'Default'};
                   end
                end
            end
            if isfield(obj.categories, 'coordinates')
                coordinates = fieldnames(obj.categories.coordinates.options);
                for k = 1:length(coordinates)
                    suggestions(end+1, :) = {'coordinates', coordinates{k}, 'Default'};
                end
                
            end
            
        end
        
        function setSolutionLabel(obj, compconf)
           obj.solutionlabel =  compconf.getSolutionLabel();
        end
        
        function l = getSolutionLabel(obj)
            l = obj.solutionlabel;
        end
    end
    
    
    
    
    
end