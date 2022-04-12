function succes = init_sbmltoolbox()
% INIT_SBMLTOOLBOX initialize SBML libraries
%
% Best is if the user installs libSBML and SBMLToolbox on the PC
% from the official sources
% see documentation
% (1) libSBML
%   http://sbml.org/Software/libSBML/docs/cpp-api/libsbml-installation.html
%   (make sure to configure the interface for MATLAB)
% (2) SBML Toolbox
%   http://sbml.org/Software/SBMLToolbox/docs/ManualV3.pdf
%   (on some systems, the installer must be run from that MATLAB command
%   line)
%
% In the case of a 32-bits windows computer, another option is to 
% download a directory win32 into matcont<version>\SBML .
% This directory contains all the necesary libaries for the import to work
% and can be obtained from http://sourceforge.net/projects/matcont/
% see files -> matcont -> matcont<version>.
%

    % try it, if it works, libsbml and SBML Toolbox is working
    succes = test_libsbml();

    if ~succes
        % no succes, but not all is over, check if libs are locally installed
        % for instance, if on windows 32 bits, arch is win32
        % so the SBML libiraries can be installed in
        % matcont<version>/SBML/win32
        % for all arch, see 'doc computer'
        
        % get /SBML dir
        sbml_dir = mfilename('fullpath');
        sbml_dir = sbml_dir(1:end-length(mfilename));
        % get /SBML/win32 or /SBML/glnx86, see 'doc computer'
        sbml_lib_dir = fullfile(sbml_dir,computer('arch'));
        if exist(sbml_lib_dir, 'dir')
            % the path with libs exists
            addpath(genpath(sbml_lib_dir));
            % now try again
            succes = test_libsbml();
        end
        % we replace Substitute.m from the SBML Toolbox, make sure
        % our Substitute is taken
        addpath(sbml_dir);
    end
    
    if ~succes
        % still no succes!
        % display some helpfull information to the user
        help('init_sbmltoolbox');
    end

end

function succes = test_libsbml()
% TEST_LIBSBML tries out a function of libsbml to look if it works

    % get sbml dir
    sbml_dir = mfilename('fullpath');
    sbml_dir = sbml_dir(1:end-length(mfilename));
    % try a libsbml function, if it works, we are ok
    succes = true;
    try
        test_file_name = fullfile(sbml_dir, 'test.xml');
        % TranslateSBML is a libSBML function
        model = TranslateSBML(test_file_name);
        % GetAllParameters is a SBML Toolbox function
        GetAllParameters(model);
    catch
        succes = false;
    end
    
end