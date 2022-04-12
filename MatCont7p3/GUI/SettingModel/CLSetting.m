classdef CLSetting < handle
   properties
       value
       validitycheck
       groupid
       subgroupid
       itemid
       displayname
       help
       visible;
       adjusted
       editable = true;
       
   end
  

   methods
       function obj = CLSetting(displayname, initvalue, checker, groupid, subgroupid, itemid, help)
           
           obj.value = initvalue;
           obj.validitycheck = checker;
           obj.groupid = groupid;
           obj.itemid = itemid;
           obj.displayname = displayname;
           
           obj.help = help;
           obj.subgroupid = subgroupid;
           obj.visible = true;
           obj.adjusted = false;
       end
       
       
       function value = getValue(obj)
           value = obj.value;
           
       end
       
       function [valid, msg] = setValue(obj, newvalue)
           if ~obj.editable
               valid = 0; msg = 'not editable';return
           end
           [valid, msg] = obj.validitycheck.validate(newvalue);
           if valid
              obj.value = newvalue; 
              obj.adjusted = true;
           end
       end
       function forceValue(obj, newvalue)
          obj.value = newvalue; 
       end
       
       function id = getGroupID(obj)
           id = obj.groupid;
       end
       function id = getSubGroupID(obj)
           id = obj.subgroupid;
       end
       function id = getItemID(obj)
           id = obj.itemid;
       end       
       
       function row = getIDs(obj)
          row = [obj.getGroupID(), obj.getSubGroupID(), obj.getItemID()]; 
       end
       
       function s = toString(obj)
          s = obj.validitycheck.toString(obj.getValue()); 
       end
       
       function newobj = copy(obj, newsettings)
           newobj = CLSetting(obj.displayname, obj.value, obj.validitycheck, obj.groupid, obj.subgroupid, obj.itemid, obj.help);
           obj.copyOver(newobj);
       end
       function h = getHelpStr(obj)
           h = obj.help; %fix?
       end
       function b = isVisible(obj)
          b = obj.visible; 
       end
       function setVisible(obj, bool)
          obj.visible = bool; 
       end
       function b = isAdjusted(obj)
           b = obj.adjusted;
       end
       
       function box = renderGUI(obj, session, settings, label, panelhandle, options, suggestions)
           GUIrender = obj.validitycheck.getGUIRenderer();
           grid = cell(1, 2);
           grid{1} = GUISettingLabel(panelhandle, settings, label, options{:});
           grid{2} = GUIrender(panelhandle, settings, label, options{:});
           
           isBool =  strcmp(obj.validitycheck.getType(), 'BOOL');
           
           box = LayoutBox(grid);
           box.fillLastAlignRight = isBool;

           box.widths(1) = max(suggestions.labelsize, box.widths(1)); 
           
       end
       
       function copyOver(obj, newobj)
            newobj.visible = obj.visible;
       end
       
       function v = getValueType(obj)
          v = obj.validitycheck.getType(); 
       end
       function n = getDisplayName(obj)
          n = obj.displayname; 
       end
       function b = sanityCheck(obj, settings)
          b = 1; 
       end
   end
    methods(Static)
       function install(settings, name, initval, restrict, catdata)
           if ~settings.activateParameter(name)
               settings.installParameter(name, CLSetting(name, initval, restrict, catdata(1), catdata(2), catdata(3), CLSettingHelp.getHelp(name) ));
           end
       end        
    end
   
    
end
