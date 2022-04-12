classdef SysGUI
    
    methods(Static)
        
        function outhandle = new()
            global path_sys
            [path_sys, ~, ~] = fileparts(which('standard.m'));
            path_sys = [path_sys filesep()];
            systems_standalone('init');
            fhandle = systems_standalone();
            if nargout == 1
               outhandle = fhandle; 
            end
        end
        
        
        function outhandle = edit(systemname)
            
            if (~ischar(systemname))
                systemname = func2str(systemname);
            end
            
            global path_sys;
            global gds;
            [path_sys, ~, ~] = fileparts(which('standard.m'));
            path_sys = [path_sys filesep()];
            load( [path_sys  systemname '.mat' ]);  %overwrites gds.
            
            fhandle = systems_standalone();
            if nargout == 1
                outhandle = fhandle;
            end
        end
        
        function outhandle = userfunctions(systemname)
            
            if (~ischar(systemname))
                systemname = func2str(systemname);
            end
            global path_sys
            global gds;
            [path_sys, ~, ~] = fileparts(which('standard.m'));
            path_sys = [path_sys '/'];
            systemmatfile = [path_sys  systemname '.mat' ];
            load(systemmatfile); %overwrites 'gds'.
            
            
            fhandle = userfun_standalone();
            if nargout == 1
                outhandle = fhandle;
            end
            
        end
        
        
        function sys = gui_loader(systemcmd)
            global gds;
            gds = [];
            fhandle = systemcmd();
       
            if isvalid(fhandle); uiwait(fhandle); end
            
            if ~isempty(gds) && isfield(gds, 'ok') && gds.ok
                sysname = gds.system;
                [path_sys, ~, ~] = fileparts(which('standard.m'));
                syspath = fullfile(path_sys, [sysname '.mat']);
                sys = CLSystem(syspath);
                
            else
                sys = [];
            end
        end
        
    end

end
