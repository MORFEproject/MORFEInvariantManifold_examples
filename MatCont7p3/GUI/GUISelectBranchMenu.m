classdef GUISelectBranchMenu < handle  

    properties
        handle
        
        eventlistener
        type;
        menuitems = []
    end
    
    methods
        function obj = GUISelectBranchMenu(parent, type,  session , branchmanager, varargin)
            if(strcmp(type,'context'))
                obj.handle = uicontextmenu;
            else
                obj.handle = uimenu(parent, 'Label' , 'Curve'  , 'Enable' , 'off' , varargin{:});
            end
            obj.type = type;
           
           obj.eventlistener = branchmanager.addlistener('newlist' , @(o,e) obj.configure(session,branchmanager));
           branchmanager.addlistener('selectionChanged' , @(o,e) obj.setSelected(branchmanager));
           
           
           obj.configure(session, branchmanager);
           obj.setSelected(branchmanager);
        end
        
        
        function configure(obj, session, branchmanager)
            hs = allchild(obj.handle);
            if (~isempty(hs)), delete(hs); end
            
            
            list = branchmanager.getCompConfigs();
            len = length(list);
            
             if (~strcmp(obj.type,'context'))
               set(obj.handle , 'Enable' , CLbool2text(len ~= 0));
             end           
            
            obj.menuitems = cell(1, len);
            
            for i = 1:length(list)
               obj.menuitems{i} = uimenu(obj.handle, 'Label' , list{i}.toString() , 'Callback' ,  @(o,e) branchmanager.select(session, i)); 
            end

        end
        function setSelected(obj, branchmanager)
           j =  branchmanager.getSelectionIndex();           
           for i = 1:length(obj.menuitems)
              set(obj.menuitems{i} , 'Checked' , CLbool2text( i == j)); 
           end
        end
    end
    
end

