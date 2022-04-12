classdef CompConf
%computation configuration    
    
   methods
       function configureSettings(~, settings)
           settings.addSetting('IP', CLSettingIP('EMPTY'));
           
           settings.addSetting('option_pause', CLSetting('Suspend Computation', 'At Special Points', InputRestrictionCategory({'At Special Points', 'At Each Point', 'Never'}), 0, 5, 1, CLSettingsHelp.getHelp('option_pause') ));
           settings.addSetting('option_archive', CLSetting('Archive Filter', 2, InputRestrictions.INT_g0, 0, 5, 2, CLSettingsHelp.getHelp('option_archive') ));
           settings.addSetting('option_output', CLSetting('Output Interval', 1, InputRestrictions.INT_g0, 0, 5, 3, CLSettingsHelp.getHelp('option_output') ));

           
           settings.addSetting('numericconfig', NumericConfig());
           settings.addSetting('previousplots', CLSettingPlotMemory());
           
       end
       
       function b = isConditional(~)
           b = false;
       end
       
       function s = getLabel(~)
           s = '<nothing>';
           
       end
       function s = toString(~)
           s = '<nothing>'; 
       end
        function s = toStringCL(obj) 
            s = obj.toString();
        end
       function s = getName(obj)
          s =  obj.toString();
       end
       function l = getSolutionLabel(obj)
          l = '  ';
       end
       
       function alist = actions(~, ~)
           alist = [];
       end
       
       function b = isAvailable(~, ~)
            b = 0;
       end
       
       function p = getPrioritynumber(~)
           p = -1;
       end
        
       function oi = getOutputInterpreter(~)
           oi = [];
       end
       
       function b = isHidden(~)
          b = 0; 
       end
       function msg = check(~, ~)
           msg = '';
       end
   end
    
end
