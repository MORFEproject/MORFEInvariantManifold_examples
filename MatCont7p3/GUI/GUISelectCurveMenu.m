classdef GUISelectCurveMenu < handle
    
    properties
        handle
        
        eventlistener
        menuitems = []
    end
    
    methods
        function obj = GUISelectCurveMenu(parent, session , branchmanager, varargin)
            
            obj.handle = uimenu(parent, 'Label' , 'Curve', 'DeleteFcn' , @(o,e) obj.destructor()  , varargin{:});
            
            %curvetypes = {'O', '-', 'EP', '-', 'LC', '-', 'H', 'LP', '-', 'PD', 'LPC', 'NS', '-', 'BP', 'BPC', '-', 'DO', '-', ...
            %    'Hom', 'HSN', '-', 'ConnecHom', 'ConnecHSN', 'HTHom', 'HTHSN', '-', 'ConnecHet', 'HTHet', 'Het'};
            curvetypes = {'O', '-', 'EP', '-', 'LC', '-', 'H', 'LP', '-', 'PD', 'LPC', 'NS', '-', 'BP', 'BPC', '-', ...
                'Hom', 'HSN', '-', 'ConnecHom', 'ConnecHSN', 'HTHom', 'HTHSN', '-', 'ConnecHet', 'HTHet', 'Het'};
            
            obj.menuitems = {};
            for k = 1:length(curvetypes)
               label = curvetypes{k}; 
               if ~strcmp(label, '-')
                   separate = k-1 >= 1 && strcmp(curvetypes{k-1}, '-');
                   fullname = DefaultValues.CURVENAMES.(label);
                   h = uimenu(obj.handle, 'Label', fullname, 'Separator', CLbool2text(separate), 'UserData', struct('label', label, 'fullname', fullname), 'Enable', 'off');
                   obj.menuitems{end+1} = h;
                   
               end
            end

            obj.eventlistener = [branchmanager.addlistener('newlist' , @(o,e) obj.configure(session,branchmanager)), ...
                                 branchmanager.addlistener('selectionChanged' , @(o,e) obj.configure(session,branchmanager))];
            obj.configure(session, branchmanager);
     
            
        end
        
        function destructor(obj)
           delete(obj.eventlistener);
           delete(obj);
            
        end
        function configure(obj, session, branchmanager)
            list = branchmanager.getCompConfigs();
            len = length(list); 
            selectedIndex =  branchmanager.getSelectionIndex();

            for i = 1:length(obj.menuitems)
                set(obj.menuitems{i}, 'Enable', 'off', 'Label', obj.menuitems{i}.UserData.fullname, 'Callback' , [], 'Checked', 'off');
                
                for branchIndex = 1:len
                    branch = list{branchIndex};
                    selected = branchIndex == selectedIndex;
                    if strcmp(obj.menuitems{i}.UserData.label, branch.getSolutionLabel()) && strcmp(obj.menuitems{i}.Enable, 'off')
                        if any(strcmp(obj.menuitems{i}.Label, {'Orbit', 'ConnectionSaddle'}))
                              set(obj.menuitems{i}, 'Enable', 'on', 'Label', branch.getName(), 'Callback' , @(o,e) branchmanager.select(session, branchIndex), 'Checked', CLbool2text(selected));
                          
                        else
                        
                            set(obj.menuitems{i}, 'Enable', 'on', 'Label', branch.toString(), 'Callback' , @(o,e) branchmanager.select(session, branchIndex), 'Checked', CLbool2text(selected));
                        end
                    end
                end
                
            end
        end
    end
    
end

