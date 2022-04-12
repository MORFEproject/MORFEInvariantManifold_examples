classdef CLSettingPlotMemory < CLSettingInterface
   
    properties
        MAXPLOTS = 8;
        plotconfigs = {};
    end
    
    methods
        
        function obj = CLSettingPlotMemory()
            obj.plotconfigs = {};
            
        end
        
        function addPlot(obj, plotconf)
            
            data = plotconf.getData();
            
            %remove plots with similar config:
            purgelist = cellfun(@(x) strcmp(x.label, data.label), obj.plotconfigs);
            obj.plotconfigs(purgelist) = [];
            
            obj.plotconfigs = [{data}; obj.plotconfigs];
            
            if length(obj.plotconfigs) > obj.MAXPLOTS
                obj.plotconfigs = obj.plotconfigs(1:obj.MAXPLOTS);
            end
            
        end
        function c = getConfigs(obj)
            c = obj.plotconfigs;
        end
        
        function plotconf = recover(obj, dimension, windowmanager)
            plotconf = [];
            for k = 1:length(obj.plotconfigs)
                if obj.plotconfigs{k}.dimension == dimension
                    plotconf = GUIPlotConf(dimension, windowmanager);
                    plotconf.restoreData(obj.plotconfigs{k});
                
                   return; 
                end
            end
        end
        
        function newobj = copy(obj, ~)
            newobj = CLSettingPlotMemory();
            newobj.plotconfigs = obj.plotconfigs;
        end        
        
    end
end