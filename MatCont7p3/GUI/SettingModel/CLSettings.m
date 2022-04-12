classdef CLSettings < handle
    properties
       fields = struct(); 
       doWarnings = 1;
       compconf;
       
       
       versionnumber = 0; %set by 'session' in order to prevent working with data from older versions of the program.

    end
    events
        settingChanged
        
    end
   methods
       function obj = CLSettings()
       end

       function bool = installSetting(obj, fieldname, paramobj)
          bool = 1;
          obj.fields.(fieldname) = paramobj;
       end
       
       function bool = addSetting(obj, fieldname, paramobj)
          bool = ~isfield(obj.fields, fieldname);
          if bool
            obj.fields.(fieldname) = paramobj;
          else
             obj.fields.(fieldname).setVisible(1); 
          end
       end       

       function disp(obj)
           
            [names, matrix] = sortedfieldnames(obj);
            maxlen = 10;
            for i = 1:length(names)
               if obj.fields.(names{i}).isVisible() && length(names{i}) > maxlen
                   maxlen = length(names{i});
               end
            end
            formatstr = '%s %s: %s\n';
    
            maincat = -1;
            subcat = 0;
            if ~isempty(names)
               subcat = matrix(1, 2);
            end
            
            hasOutput = false;
            hasOutput2 = false;
            for i = 1:length(names)
               item = obj.fields.(names{i});
               if maincat ~= matrix(i, 1)
                  maincat = matrix(i, 1);
                  subcat = matrix(i, 2);
                  if hasOutput
                    fprintf('%s\n', repmat('=', 1, 1*maxlen + 10));
                    hasOutput=false;
                  end
               end
               if subcat ~= matrix(i, 2)
                  subcat = matrix(i, 2);
                  if hasOutput2
                    fprintf('\n');
                    hasOutput2 = false;
                  end
               end
               if item.isVisible()
                   namehtml = ['<a href="matlab: disp(''' item.getHelpStr() ''')">' names{i} '</a>'];
                   fprintf(formatstr, repmat(' ', 1, maxlen - length(names{i})), namehtml , item.toString());
                   hasOutput = true;
                   hasOutput2 = true;
               end
            end

       end
       
       function catlist = getSettings(obj)
           [names, matrix] = sortedfieldnames(obj);
           items = cell(1, length(names));
           for i = 1:length(names)
               items{i} = obj.fields.(names{i});
               
           end
           mask = find(cellfun(@(x) 0+x.isVisible(), items));
           names = names(mask);
           matrix = matrix(mask, :);
           catlist = {};
           for maincat = unique(matrix(:,1))'
                mask = find(matrix(:,1) == maincat);
                subcatlist = {};
                for subcat = unique(matrix(mask,2))'
                    subcatlist{end+1} = [num2str(subcat);names(mask(matrix(mask,2) == subcat))];
                end
                catlist{end+1} = [maincat, subcatlist];
           end

           
           
       end
       
       function [bool, warnmsg] = setValue(obj, fieldname, value)
           if obj.exists(fieldname)
                [bool, warnmsg] = obj.fields.(fieldname).setValue(value);
                
                if ~bool && obj.doWarnings
                   fprintf(2, ['ERROR ' fieldname ': ' warnmsg '\n']); 
                elseif bool
                   obj.notify('settingChanged'); 
                end
           end
       end

       function value = getValue(obj, fieldname)
           if obj.exists(fieldname)
               value = obj.fields.(fieldname).getValue();
           end
       end
       function value = getVal(obj, fieldname, defaultvalue)
           if isfield(obj.fields, fieldname) && obj.fields.(fieldname).isVisible()
               value = obj.fields.(fieldname).getValue();
           else
               value = defaultvalue;
           end
       end       
       function bool = exists(obj, fieldname)
           bool = isfield(obj.fields, fieldname);
           if ~bool && obj.doWarnings
              fprintf(2, ['"' fieldname '"' ' does not exist!\n']);
           end
       end
      
       function s = getSetting(obj, fieldname)
           if isfield(obj.fields, fieldname)
               s = obj.fields.(fieldname);
           else
              s = []; 
           end
       end
       function s = setSetting(obj, fieldname, setting)
           obj.fields.(fieldname) = setting;
       end
       
       function b = refresh(obj)
          b = 1; obj.notify('settingChanged'); 
       end
       
       function out = subsref(obj , S)               
               if isfield(obj.fields, S(1).subs)
                  %disp('intercept');
                  if length(S) == 1
                      out = obj.fields.(S(1).subs).getValue();
                  
                  elseif length(S) == 3 && strcmp(S(3).type, '()') && strcmp(S(2).subs, 'set') && length(S(3).subs) == 1
                        value = S(3).subs{1};
                        setValue(obj, S(1).subs, value);
                  else
                      fprintf(2, 'ERROR\n');
                  end
               else
               
                    out = builtin('subsref', obj, S);
               end
       end
       
       function newobj = copy(obj)
           settingnames = obj.sortedfieldnames();
           newobj = CLSettings();
           
           for index = 1:length(settingnames)
              name = settingnames{index}; 
              setting = obj.fields.(name);
              copysetting = setting.copy(newobj);
              if ~isempty(copysetting)
                 newobj.setSetting(name, copysetting); 
              end
              
           end
           
       end
       function b = setVisible(obj, settingname, bool)
            obj.fields.(settingname).setVisible(bool); b = 1;
       end
       
       function box = renderGUI(obj, session, fieldname, panelhandle, options, suggestions)
            box = obj.fields.(fieldname).renderGUI(session, obj, fieldname, panelhandle, options, suggestions);
       end
       function o = gui(obj, index)
           if nargin < 2
              index = 1; 
           end
          o =  GUIRenderSettings(obj, index);
       end
       
       function b = sanityCheck(obj)
          names = fieldnames(obj.fields);
          b = 1;
          for k = 1:length(names)
              if isfield(obj.fields, names{k}) && ~obj.fields.(names{k}).sanityCheck(obj)
                    obj.removeSetting(names{k});
                    b = 0;
              end
          end
          
          if ~b
             fprintf(2, 'WARNING: The system seems to have changed (coord, param or userfunction mismatch), some settings were resetted\n'); 
              
          end
          
           
       end
       function b = removeSetting(obj, settingname)
           b = 0;
           if isfield(obj.fields, settingname)
               %fprintf('CLSETTINGS: removing %s\n', settingname);
               delete(obj.fields.(settingname));
               obj.fields = rmfield(obj.fields, settingname);
               b = 1;
           end
       end
       
   end

   methods(Hidden)
       function [sortednames, matrix] = sortedfieldnames(obj)
           
           names = fieldnames(obj.fields);
           matrix = zeros(length(names), 3);
           settings = cell(length(names), 1);
           for i = 1:length(names)
               setting = obj.fields.(names{i});
               settings{i} = setting;
               matrix(i,:) = setting.getIDs();

           end
           [matrix, index] = sortrows(matrix);
           sortednames = names(index);
           %{
           for i = 1:length(sortednames)
              fprintf('%s  %s\n', sortednames{i},  mat2str(matrix(i, :))); 
               
           end
           %}
       end
       function b = test(obj)  %FIXME TODO REMOVE
           b = 1;
           names = fieldnames(obj.fields);
           matrix = zeros(length(names), 3);
           settings = cell(length(names), 1);
           for i = 1:length(names)
               setting = obj.fields.(names{i});
               settings{i} = setting;
               matrix(i,:) = setting.getIDs();
               
           end
           [matrix, index] = sortrows(matrix);
           sortednames = names(index);
           for i = 1:length(sortednames)
               label = sprintf('s%s_%s', num2str(matrix(i, 1)), num2str(matrix(i, 2)));
               if isfield(DefaultValues.SECTIONNAMES, label)
                    label = DefaultValues.SECTIONNAMES.(label);
               end
               
               fprintf('%s  %s  %s\n', pad(sortednames{i}, 20),  pad(mat2str(matrix(i, :)), 20), label)
           end
       end
       
       function b = testHelp(obj) %FIXME TODO REMOVE
           names = fieldnames(obj.fields);
           for x = 1:length(names)
              if isempty(CLSettingsHelp.getHelp(names{x})) 
               fprintf('%s needs help\n', names{x});
              end
           end
           b = 1;
           
       end
       
       
       function out = softClear(obj)
           names = fieldnames(obj.fields);
           for i = 1:length(names)
               obj.fields.(names{i}).setVisible(0);
           end
           out = [];
       end
       function out = revealAll(obj)
           names = fieldnames(obj.fields);
           for i = 1:length(names)
               obj.fields.(names{i}).setVisible(1);
           end
           out = [];
       end       
       
   end

end   
  
