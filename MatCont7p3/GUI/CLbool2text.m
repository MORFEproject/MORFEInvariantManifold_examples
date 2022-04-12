function s = CLbool2text(bool)
% true/false -> on/off, for use in MATLAB GUI. %.

    if bool
       s = 'on'; 
    else
       s = 'off';
    end
end
