function init()


%only disable these two warnings during continuation
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');

% turn off all warnings: (old setting)
%warning off;



[list,val] = spparms;
spparms('default');
for i=1:length(list(:,1))
    if strcmp(list(i,:),'umfpack')
%             spparms('umfpack',0); %switch umfpack off, use v4solver
    end
end

addpath(cd);
addpath([cd '/Continuer/']);
addpath([cd '/Equilibrium/']);
addpath([cd '/LimitCycle/']);
addpath([cd '/PeriodDoubling/']);
addpath([cd '/Systems/']);
addpath([cd '/LimitPoint/']);
addpath([cd '/Hopf/']);
addpath([cd '/LimitPointCycle/']);
addpath([cd '/NeimarkSacker/']);
addpath([cd '/BranchPoint/']);
addpath([cd '/BranchPointCycle/']);
addpath([cd '/Homoclinic/']);
addpath([cd '/HomoclinicSaddleNode/']);
addpath([cd '/HomotopySaddle/']);
addpath([cd '/HomotopySaddleNode/']);
addpath([cd '/HomotopyHet/']);
addpath([cd '/Heteroclinic/']);
addpath([cd '/MultilinearForms/']);
addpath([cd '/Help/']);
addpath([cd '/LimitCycleCodim2/']);
addpath([cd '/SBML/']);
addpath([cd '/Testruns/']);
addpath([cd '/Testruns/TestSystems/']);

%    addpath([cd '/Help/']);
% Find and go to correct directory (class directory)
p = mfilename('fullpath');
p = p(1:length(p)-length(mfilename));
p = strcat(p,'/LimitCycle');
curdir = cd;
cd(p);

% Compile the c-files (optimized)
if ~(exist(strcat('BVP_LC_jac.',mexext),'file'))
if ~isempty(regexp(mexext,'64','match'))
    mex -largeArrayDims -O BVP_LC_jac.c;
    mex -largeArrayDims -O BVP_PD_jac.c;
    mex -largeArrayDims -O BVP_BPC_jacC.c;
    mex -largeArrayDims -O BVP_BPC_jacCC.c;
    mex -largeArrayDims -O BVP_LPC_jac.c;
    mex -largeArrayDims -O BVP_NS_jac.c;
    mex -largeArrayDims -O BVP_LCX_jac.c;
else
    mex -O BVP_LC_jac.c;
    mex -O BVP_PD_jac.c;
    mex -O BVP_BPC_jacC.c;
    mex -O BVP_BPC_jacCC.c;
    mex -O BVP_LPC_jac.c;
    mex -O BVP_NS_jac.c;
    mex -O BVP_LCX_jac.c;
end
end


% Return to directory we started in
cd (curdir);

global MC;
MC = [];
clear MC;
