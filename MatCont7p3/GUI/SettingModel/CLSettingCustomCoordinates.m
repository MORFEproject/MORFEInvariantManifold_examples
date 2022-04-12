classdef CLSettingCustomCoordinates < CLSetting

    properties
       coordinates 
       prefix
       tag
    end
    
    methods
        function obj = CLSettingCustomCoordinates(settings, coordinates, tag, prefix, subcategory)
            dim = length(coordinates);
            obj = obj@CLSetting(tag, zeros(1, dim), InputRestrictions.vector(dim), 2, subcategory, length(coordinates) + 1, '~~~');
            obj.coordinates = coordinates;      
            obj.prefix = prefix;
            if ~isempty(settings)
                obj.installGhosts(settings);
            end
            obj.subgroupid = subcategory;
            obj.tag = tag;
        end
        function box = renderGUI(varargin)
            box = [];
            
        end
        
        function installGhosts(obj, settings)
            for i = 1:length(obj.coordinates)
                coordinate = obj.coordinates{i};
                settings.addSetting([obj.prefix, coordinate], CLSettingCoordinate(obj, i, obj.subgroupid));
            end
            
            
        end
        
       function newobj = copy(obj, newsettings)
           newobj = CLSettingCustomCoordinates(newsettings, obj.coordinates, obj.tag, obj.prefix, obj.subgroupid);
           newobj.value = obj.value;
           obj.copyOver(newobj);
       end        
        
       function b = sanityCheck(obj, settings)
            b = obj.internalSanityCheck(settings);
           
            if ~b
                for i = 1:length(obj.coordinates)
                    coordinate = obj.coordinates{i};
                    settings.removeSetting([obj.prefix, coordinate]);
                end
            end
            
       end
       
       function b = internalSanityCheck(obj, settings)
          system = settings.system;
          if isempty(system); b = 0; return; end
          coordinates = system.getCoordinates();
          if length(coordinates) ~= length(obj.coordinates); b = 0; return; end
          
          for k = 1:length(obj.coordinates)
                if ~strcmp(obj.coordinates{k}, coordinates{k})
                    b = 0;
                    return
                end
          end
          b = 1;
       end          
    end
end
