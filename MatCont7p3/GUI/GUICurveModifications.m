classdef GUICurveModifications
    %Adds special plot-options to certain curves under certain conditions
    %for continuer only.
    %The Continuer-Outputter (GUIPlotOutputter) will use this class to
    %determine if any modifiers are avaible for a plot. 
    % A modifier-function retrieves continuer data and can choose to add
    % plot-options to the current curve or to leave the curve unmodified
    % (by returning {})
    % For every modifier, GUILineConfigPanel will add a field for the user
    % to fill in. A non-empty field activates the modifier.
    
    
    properties(Constant)
        
        %structure with the modification. Add 
        MODS = struct(...
            'EP', struct('title', 'if unstable', 'function', @EPmod),... %Equilibrium
            'LC', struct('title', 'if unstable', 'function', @LCmod));   %Limit Cycle
        
        
    end

    methods(Static)
    
        function mod = getMod(curvelabel, lineoptions)
            %returns a modifier function if a modifier exist and if the
            %user has entered a plot-option in the GUILineConfigPanel.
            
            mod = {};
            
            %obtain plot-option from the line options ('Plot
            %Properties'/GUILineConfigPanel). Entered by user.
            optdata = lineoptions.getLineModSetting(curvelabel);
            
            %user entered nothing: modifier ignored
            if isempty(optdata); return; end
            
            %mod exists, user has entered plot-opts, return a mod-function
            if GUICurveModifications.hasMod(curvelabel)
               data =  GUICurveModifications.MODS.(curvelabel);
               mod = @(contdata, index) data.function(contdata, index, optdata);
            end
            
        end
        
        %Check if modifier has been configured, return name.
        function b = hasMod(curvelabel)
            if isfield(GUICurveModifications.MODS, curvelabel)
               data =  GUICurveModifications.MODS.(curvelabel);
               b = data.title;
            else
               b = []; 
            end
        end
        
        
    end
end



function mod = EPmod(contdata, index, plotopts)
%Modifier for Equilibrium. Returns 'plotopts' (entered by user) if
%unstable.
%contdata: {x,h,f}
f = contdata{3};
%f has non-Nan values and not all eigenvalues have real-part smaller than 0
if ~any(isnan(f(:,index))) && ~all(real(f(:,index)) < 0)  
    mod = plotopts; %unstable: add  the settings
else
    mod = {}; %stable: do nothing
end
end

function mod = LCmod(contdata, index, plotopts)
%Modifier for Limit Cycle. Returns 'plotopts' (entered by user) if
%unstable.
%contdata: {x,h,f}
f = contdata{3};
%skip mesh '1' to skip mesh.
s = find(abs(f(:, index) - 1) < 1e-12, 1) + 1;

mult = f(s:end, index); %extract multipliers
[~, index] = min(abs(mult - 1)); %find multiplier 1 and remove
mult(index) = [];

if any(isnan(f(:,index))) || all(abs(mult) < 1)
    mod = {};  %if stable or disabled, don't modify settings
else
    mod = plotopts;  %if unstable, add new settings
end
end