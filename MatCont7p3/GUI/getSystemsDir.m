function d = getSystemsDir()
d = fullfile(fileparts(fileparts(which('MATCONTGUI.m'))), 'Systems');

end