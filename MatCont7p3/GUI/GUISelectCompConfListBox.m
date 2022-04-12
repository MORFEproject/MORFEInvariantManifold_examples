classdef GUISelectCompConfListBox < handle  

    
    properties
        handle
        eventlistener

        

        labels
        listener1
    end
    
    methods
        function obj = GUISelectCompConfListBox(parent, session, branchmanager, names, labels, varargin)

           
            obj.labels = labels;
            
            
            obj.handle = uicontrol(parent, 'Style' , 'popupmenu', 'String', names,   varargin{:});
            obj.listener1 = branchmanager.addlistener('selectionChanged' , @(o,e) obj.setSelected(branchmanager));
            set(obj.handle , 'Callback' , @(o,e) obj.callback(session , branchmanager), 'DeleteFcn', @(o, e) obj.destructor());
            
            obj.setSelected(branchmanager);
        end
        
        function callback(obj, session, branchmanager)
            index = get(obj.handle, 'Value');
            
            branchmanager.setConf(session, obj.labels{index});
            
        end
        
        function setSelected(obj, branchmanager)
            selectlabel = branchmanager.getSelection().getLabel();
            
            index = find(strcmp(obj.labels, selectlabel), 1);
            
            if isempty(index) %mismatch with branchmanager, disabling..
                set(obj.handle, 'String', {'   '}, 'Value', 1, 'Enable', 'off');
                obj.labels = {};
            else
                set(obj.handle , 'Value' , index); 
                
            end
            
        end
        
        function destructor(obj)
           delete(obj.listener1);
           delete(obj);
        end
        
        
        
    end
    
end
