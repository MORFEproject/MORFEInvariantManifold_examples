function matcont(varargin)
% Start up Matcont-GUI
%
% >>> matcont
% no arguments: GUI will load previous state
%
% >>> matcont clean
% one argument: GUI will load an empty state
%
% >>> matcont clean reset
% two arguments: GUI will first 'clear all' and then load an empty state to avoid any cached code.


MINIMALVERSION = '9.2'; 
VERSIONNAME = '2017a';

if (verLessThan('matlab', MINIMALVERSION))
    for i = 1:1
        fprintf(2, 'matlab version needs to be %s (%s) or higher\n', MINIMALVERSION, VERSIONNAME);
    end
    pause(2);
end




%Same init as CL-version.
init();
addpath('GUI');


workingpath = pwd(); %save working dir

try
%execute path-init of GUI
cd(fullfile(workingpath, 'GUI'));
initpath();

cd(workingpath); %restore working dir
catch
    %failure is an option if the script was run succesfully before, GUI can
    %then be started from anywhere. not just main-directory.
end 

%Start GUI
MATCONTGUI(nargin);

end
