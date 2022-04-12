function d = correctSystemsDir(path)
    d = path;
    splitted = strsplit(path, '\');
    i = find(strcmp(splitted, 'Systems'));
    
    if ~isempty(i)
        i = i(end);
        d = fullfile( getSystemsDir() , strjoin(splitted(i+1:end), filesep));
    
    else
        splitted = strsplit(path, '/');
        i = find(strcmp(splitted, 'Systems'));
        if ~isempty(i)
            i = i(end);
            d = fullfile( getSystemsDir() , strjoin(splitted(i+1:end), filesep));
            
        end
    end
    
         

end