classdef GUISelectBranchListbox < handle  

    
    properties
        handle
        eventlistener
        strlist;
        hidelist;
        displayfunction;
    end
    
    methods
        function obj = GUISelectBranchListbox(parent, session, displayfunction, varargin)
            
            branchmanager = session.branchmanager;
            obj.handle = uicontrol(parent, 'Style' , 'popupmenu' ,   varargin{:});
           obj.displayfunction = displayfunction;
            
           branchmanager.addlistener('newlist' , @(o,e) obj.configure(session,branchmanager));
           branchmanager.addlistener('selectionChanged' , @(o,e) obj.setSelected(branchmanager));
           
           obj.configure(session, branchmanager);
           
           set(obj.handle , 'Callback' , @(o,e) obj.callback(o,session , branchmanager));
        end
        
        
        function configure(obj, session, branchmanager)
            
            list = branchmanager.getCompConfigs();
           
            len = length(list);
            
            set(obj.handle , 'Enable' , CLbool2text(len ~= 0));
           
            obj.strlist = cell(1, len);
            obj.hidelist = cell(1, len);
            for i = 1:len
               obj.strlist{i} = obj.displayfunction(list{i}); 
               obj.hidelist{i} = list{i}.isHidden();
            end
            
            obj.renderList(branchmanager.getSelectionIndex(), branchmanager);
            
        end
        function renderList(obj, index, branchmanager)
            if isempty(obj.strlist)
                set(obj.handle , 'Enable' , 'off' , 'String' , {' '} , 'Value' , 1);
            else
                displaylist = {};
                indexlist = [];
                newindex = [];
                for i = 1:length(obj.strlist)
                    if ~obj.hidelist{i} || i == index
                        displaylist{end+1} = obj.strlist{i};
                        indexlist(end+1) = i;
                        if i == index
                           newindex = length(displaylist); 
                        end
                    end
                end
                if (branchmanager.validSelectionMade())
                    set(obj.handle , 'Enable' , 'on' , 'String' , displaylist, 'UserData', indexlist, 'Value', newindex);
                else
                    set(obj.handle , 'Enable' , 'on' , 'String' , [displaylist, {'     '}], 'UserData', indexlist, 'Value', length(displaylist)+1);
                end
            end
            
            
        end
        
        
        function setSelected(obj, branchmanager)
           obj.renderList(branchmanager.getSelectionIndex(), branchmanager);
        end
        
        
        function callback(obj , handle ,session, branchmanager)
            j =  get(handle, 'Value');
            map = get(handle, 'UserData');
            branchmanager.select(session, map(j));
        end
        
    end
    
end
