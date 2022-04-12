classdef GUIPlotConfLite
   properties
        dimension
        configs
        label
        figureposition;
   end
    
   methods
       function obj = GUIPlotConfLite(dimension, configs, label, position)
           obj.dimension = dimension;
           obj.configs = configs;
           obj.label = label;
           obj.figureposition = position;
           
       end
       
       
       function plotconf = reconstruct(obj)
       end
       
       function l = toString(obj)
            l = [obj.label ' - '];
            for k = 1:obj.dimension
                l = strcat(l, ['['  strip(sprintf('%.4g ', obj.configs{k}.region))   '],']);
            end
            l = strip(l, ',');
       end
   end
end