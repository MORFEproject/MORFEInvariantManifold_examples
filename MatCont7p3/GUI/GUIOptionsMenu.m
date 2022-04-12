classdef GUIOptionsMenu

    
    properties
    end
    
    methods(Static)
        function install(parent, session, optionnames)
            settings = CLSettings();
            cc = CompConf();
            cc.configureSettings(settings);
            
            for k = 1:length(optionnames)
                optionname = optionnames{k};
                setting = settings.getSetting(optionname);
                type = setting.getValueType();
                
                switch(type)
                    case 'VALUE'
                        uimenu(parent , 'Label' , setting.getDisplayName() , 'Callback' , @(o,e) valueCallback(o, session, optionname));
                    case 'BOOL'
                        uimenu(parent , 'Label' , setting.getDisplayName() , 'Callback' , @(o,e) boollistCallback(o, session , optionname)); 
                    case 'CATEGORY'
                        uimenu(parent , 'Label' , setting.getDisplayName() , 'Callback' , @(o,e) listCallback(o, session , optionname)); 
                    otherwise
                        fprintf(2, 'ERROR %s: %s\n', mfilename, optionname);
                end
            end
            
            delete(settings);
        end
        
    end
    
end




function valueCallback(handle , session , optionname)
    setting = session.settings.getSetting(optionname);
    subtitle = setting.getHelpStr();
    answer = inputdlg(subtitle, handle.Label ,1 ,{ setting.toString() });
    
    
    if (~isempty(answer))
        [valid, msg] = setting.setValue(evalin('base', answer{1}));
        
        if ~valid
           %fprintf(2, '%s: %s\n', handle.Label, msg);
           
           valueCallback(handle, session , optionname); 
        end
    end

end


function listCallback(handle , session , optionname)
    setting = session.settings.getSetting(optionname);
    subtitle = setting.getHelpStr();
    list = setting.validitycheck.categories;
    initvalue = find(strcmp(setting.getValue(), list));
    
    [selection,ok] = listdlg('ListString', list ,'SelectionMode' , 'single' , ...
        'InitialValue' , initvalue  , 'Name' , handle.Label, 'PromptString' , subtitle ,'ListSize' , [300 100] );
    
    if (ok)
       setting.setValue(list{selection}); 
    end
end

function boollistCallback(handle , session , optionname)
    setting = session.settings.getSetting(optionname);
    subtitle = setting.getHelpStr();
    list = {'false', 'true'};
    initvalue = setting.getValue() + 1; %MATLAB counts from 1: so false<->1 and true<->2
    
    [selection,ok] = listdlg('ListString', list ,'SelectionMode' , 'single' , ...
        'InitialValue' , initvalue  , 'Name' , handle.Label, 'PromptString' , subtitle ,'ListSize' , [300 100] );
    
    if (ok)
       setting.setValue(selection==2); 
    end
end


