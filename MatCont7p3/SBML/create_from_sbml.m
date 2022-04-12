function [system_m_file, gds_struct] = create_from_sbml(model)
%CREATE_FROM_SBML Using a SBML_model structure, create file for matcont
%
% this is replaced by gds_from_sbml, but I leave this code in, as it
% creates the file <system_name>.m, useable for ODE solving. In the time
% that it was written, some stuff has been added to gds_from_sbml.m that
% might be usefull here to.
%
% If this is ever needed to have a commandline only tool to make
% <system_name>.m files. Best is to use gds_from_sbml.m to get the gds
% structure, and then use that gds structure to write the file.
%
% This file can probably be safely deleted, and was written as a prototype,
% it is not fully replaced by gds_from_sbml.m

if (~isSBML_Model(model))
    error('create_from_sbml(): argument must be an SBML Model');
end

if (~isempty(model.id))
    model_name = model.id;
elseif (~isempty(model.name))
    model_name = model.name;
else
    error('create_from_sbml(): SBML Model must have a name');
end

[parameters_names, parameters_values] = GetAllParameters(model);
species = AnalyseSpecies(model);
n_species = length(species);
reactions = model.reaction;
compartments = model.compartment;
if (~isempty(model.time_symbol))
    time_symbol = model.time_symbol;
else
    time_symbol = 't';
end

% MATCONT uses coordinates, SBML uses species
coordinates_names = [];
coordinates_string = '[';
coordinates_equations = [];
coordinates_initial_values = [];
assignments_names = [];
assignments_formulas = [];
for i = 1:n_species
    specie_name = species(i).Name;
    if species(i).ChangedByReaction
        coordinates_names = [coordinates_names specie_name];
        coordinates_equations = [coordinates_equations ...
            species(i).KineticLaw];
        coordinates_initial_values = [coordinates_initial_values ...
            species(i).initialValue];
        coordinates_string = strcat(coordinates_string,specie_name,',');
    elseif species(i).ChangedByAssignmentRule
        assignments_names = [assignments_names specie_name];
        assignments_formulas = [assignments_formulas ...
            species(i).AssignmentRule];
    else % something else, handle it like a parameter
        parameters_names = [parameters_names specie_name];
        parameters_values = [parameters_values species(i).initialValue];
    end
    syms(char(specie_name)); % enter coordinate as symbol
end
coordinates_string = strcat(coordinates_string,']');
n_coordinates = length(coordinates_names);
assert(n_coordinates > 0, 'create_from_sbml(): Model has to have coordinates');
% build a matlab symbolic string of the parameters like: '[p_1,p_2,...]'
n_parameters = length(parameters_names);
parameter_string = '[';
for parameter = parameters_names
    parameter_string = strcat(parameter_string,parameter,',');
    syms(char(parameter)); % enter parameter as symbol
end
parameter_string = strcat(parameter_string,']');
% ignore the compartments, they have a 'size' set, important for matcont???
for compartment = compartments
    compartment_name = compartment.id;
    coordinates_equations = strrep(coordinates_equations,...
        compartment_name,'1.0');
end
% sbml adds the name of the reaction to the
% prameter if there might be overlap, but in
% matcont this overlap is needed
for reaction = reactions
    react_par_names = GetParameterFromReaction(reaction);
    react_unique_par_names = GetParameterFromReactionUnique(reaction);
    for i = 1:length(react_par_names)
        coordinates_equations = strrep(coordinates_equations, ...
            char(react_unique_par_names(i)), ...
            char(react_par_names(i)));
    end
end
% put 'dydt = [ .... ]' as a symbol in matlab
dydt_string = 'dydt=[ ';
for i = 1:n_coordinates
    dydt_string = strcat(dydt_string,char(coordinates_equations(i)),';');
end
dydt_string = strcat(dydt_string,'];');
eval(dydt_string); % enter dydt in the workspace

% build gds
gds_struct = default_gds();
gds_struct.system = model_name;
for i = 1:n_coordinates
    gds_struct.coordinates{i,1} = char(coordinates_names(i));
    gds_struct.coordinates{i,2} = coordinates_initial_values(i);
end
for i = 1:n_parameters
    gds_struct.parameters{i,1} = char(parameters_names(i));
    gds_struct.parameters{i,2} = parameters_values(i);
end
gds_struct.time = {char(time_symbol),0};
for i = 1:n_coordinates
    lhs = [char(coordinates_names(i)) ''''];
    rhs = char(eval(char(coordinates_equations(i))));
    rhs = substitute_assignments(rhs,...
        assignments_names,assignments_formulas);
    gds_struct.equations = strvcat(gds_struct.equations,[lhs '=' rhs]);
end
gds_struct.dim = n_coordinates;

% ------------------------------------------------------------------------
% Create file with name of the model apended with .m
system_m_file = fullfile([model_name, '.m']);
fileID = fopen(system_m_file, 'w');

% ------------------------------------------------------------------------
% Write boilerplate code
fprintf(fileID, 'function out = %s\n', model_name);
fprintf(fileID, '\n');
fprintf(fileID, 'out{1} = @init;\n');
fprintf(fileID, 'out{2} = @fun_eval;\n');
fprintf(fileID, 'out{3} = @jacobian;\n');
fprintf(fileID, 'out{4} = @jacobianp;\n');
fprintf(fileID, 'out{5} = @hessians;\n');
fprintf(fileID, 'out{6} = @hessiansp;\n');
fprintf(fileID, 'out{7} = @der3;\n');
fprintf(fileID, 'out{8} = [];\n');
fprintf(fileID, 'out{9} = [];\n');

% ------------------------------------------------------------------------
% Write the evaluation function
fprint_functionheader(fileID,'fun_eval','dydt', ...
                      assignments_names,assignments_formulas,...
                      coordinates_names,parameters_names);
fprintf(fileID, ['dydt = ' symmat2str(dydt) ';\n']);

% ------------------------------------------------------------------------
% Write the init function
fprintf(fileID, '\n');
fprintf(fileID, ['%% ----------------------------------------------' ...
                 '----------------------------\n']);
fprintf(fileID, 'function [tspan,y0,options] = init\n');
fprintf(fileID, 'handles = feval(%s)\n', model_name);
% ??? seems good to store the initial values in y0 ???
fprintf(fileID, 'y0=[');
for coordinates_value = coordinates_initial_values
    fprintf(fileID, '%g,',coordinates_value);
end
fprintf(fileID, '];\n');
fprintf(fileID, ['options = odeset(''Jacobian'',handles(3),' ...
                                  '''JacobianP'',handles(4), ...\n' ...
                                  '\t\t\t\t''Hessians'',handles(5),' ...
                                  '''HessiansP'',handles(6));\n']);
fprintf(fileID, 'tspan = [0 10];\n');

% ------------------------------------------------------------------------
% Write the jacobian function
fprint_functionheader(fileID,'jacobian','jac', ...
                      assignments_names,assignments_formulas,...
                      coordinates_names,parameters_names);
jac = jacobian(dydt,eval(char(coordinates_string)));
fprintf(fileID, ['jac = ' symmat2str(jac) ';\n']);

% ------------------------------------------------------------------------
% Write the jacobianp function
fprint_functionheader(fileID,'jacobianp','jacp', ...
                      assignments_names,assignments_formulas,...
                      coordinates_names,parameters_names);
jacp = jacobian(dydt,eval(char(parameter_string)));
fprintf(fileID, ['jac = ' symmat2str(jacp) ';\n']);

% ------------------------------------------------------------------------
% Write the hessians
fprint_functionheader(fileID,'hessians','hess',...
                        assignments_names,assignments_formulas,...
                        coordinates_names,parameters_names);
for i = 1:n_coordinates
    coordinate = coordinates_names(i);
    hess = diff(jac,eval(char(coordinate)));
    hess_string = sprintf('hess(:,:,%d) = ',i);
    fprintf(fileID,[hess_string symmat2str(hess) ';\n']);
end

% ------------------------------------------------------------------------
% Write the hessiansp
fprint_functionheader(fileID,'hessiansp','hessp',...
                        assignments_names,assignments_formulas,...
                        coordinates_names,parameters_names);
for i = 1:n_parameters
    parameter = parameters_names(i);
    hessp = diff(jac,eval(char(parameter)));
    hessp_string = sprintf('hessp(:,:,%d) = ',i);
    fprintf(fileID,[hessp_string symmat2str(hessp) ';\n']);
end

% ------------------------------------------------------------------------
% Write the der3
fprint_functionheader(fileID,'der3','tens',...
                        assignments_names,assignments_formulas,...
                        coordinates_names,parameters_names);
for i = 1:n_coordinates
    first_coordinate = coordinates_names(i);
    first_hess = diff(jac,eval(char(first_coordinate)));
    for j = 1:n_coordinates
        second_coordinate = coordinates_names(j);
        tens = diff(first_hess, eval(char(second_coordinate)));
        tens_string = sprintf('tens(:,:,%d,%d) = ',i,j);
        fprintf(fileID,[tens_string symmat2str(tens) ';\n']);
    end
end

% ??? can we do something with the parameter values?
%fprintf(fileID, '%% are the parameter values usable?\n');
%for i = 1:n_parameters
%    fprintf(fileID, '%% %s = %g;\n', char(parameters_names(i)), ...
%                                     parameters_values(i));
%end

% ------------------------------------------------------------------------
% Cleanup
fclose(fileID);

end

function fprint_functionheader(fileID, function_string, return_string, ...
                               assignments_names,assignments_formulas,...
                               coordinate_names,parameters_names)
% FPRINT_FUNCTIONHEADER
% prints to fileID:
% %% -----------------
% function <return_string> = <function_string>(t,kmrgd, par1, par2, ...)
% 
fprintf(fileID, '\n');
fprintf(fileID, ['%% -------------------------------------------' ...
                 '-------------------------------\n']);
fprintf(fileID, 'function %s = %s(t,kmrgd',return_string,function_string);
for parameter_name = parameters_names
    fprintf(fileID, ',%s', char(parameter_name));
end
fprintf(fileID, ')\n');
% matcont gives array of values, convert coordinates to the array
for i = 1:length(coordinate_names)
    fprintf(fileID, '%s = kmrgd(%d);\n', char(coordinate_names(i)), i);
end
% let MATLAB handle assignments
for i = 1:length(assignments_names)
    lhs = char(assignments_names(i));
    rhs = char(assignments_formulas(i));
    fprintf(fileID, '%s = %s;\n', lhs, rhs);
end
end

function s = symmat2str(symmat)
    [m,n] = size(symmat);
    s = '[ ...\n\t';
    for i = 1:m
        for j = 1:n
            s = strcat(s,char(symmat(i,j)),',');
        end
        s = strcat(s,'; ...\n\t');
    end
    s = strcat(s,']');
end

function formula = substitute_assignments(formula,...
    assignments_names,assignments_formulas)
% SUBSTITUTE_ASSIGNMENTS takes a formula and keeps substituting the
%                        assignments till impossible
rule_applied = 1;
iterations_left = length(assignments_names) + 1;
while rule_applied > 0 && iterations_left > 0
    rule_applied = 0;
    for i = 1:length(assignments_names)
        str = formula;
        exp = ['\<' assignments_names(i) '\>'];
        repstr = assignments_formulas(i);
        formula = regexprep(str,exp,repstr);
        rule_applied = rule_applied + strcmp(str, formula)==false;
    end
    iterations_left = iterations_left - 1;
end
assert(rule_applied == 0, ...
    'substitute_assignments(): Cyclic dependency of rules dedected');
end

function gds = default_gds()
% DEFAULT_GDS returns an empty gds structure
% TODO: the gui code in GUI/systems.m could make use of this
    gds = [];
    gds.coordinates = [];
    gds.parameters = [];
    gds.time = {'t' 0};
    gds.options.InitStepsize = [];
    gds.options.MinStepsize = [];
    gds.options.MaxStepsize = [];
    gds.options.MaxCorrIters = [];
    gds.options.MaxNewtonIters = [];
    gds.options.MaxTestIters = [];
    gds.options.MoorePenrose = [];
    gds.options.SymDerivative = [];
    gds.options.SymDerivativeP = [];
    gds.options.Increment = [];
    gds.options.FunTolerance = [];
    gds.options.VarTolerance = [];
    gds.options.TestTolerance = [];
    gds.options.Singularities = [];
    gds.options.MaxNumPoints = [];
    gds.options.Backward = [];
    gds.options.CheckClosed = [];
    gds.options.TestFunctions = [];
    gds.options.WorkSpace = [];
    gds.options.Locators = [];
    gds.options.Adapt = [];
    gds.options.IgnoreSingularity = [];
    gds.options.ActiveParams = [];
    gds.options.Multipliers = [];
    gds.options.Eigenvalues = [];
    gds.options.Userfunctions = [];
    gds.options.UserfunctionsInfo = [];
    gds.options.PRC = 0;
    gds.options.dPRC = 0;
    gds.options.Input = [];
    gds.options.ActiveUParams = [];
    gds.options.ActiveSParams = [];
    gds.options.ActiveSParam = [];
    gds.options.ActiveParams = [];
    gds.system = '';
    gds.curve.new = '';
    gds.curve.old = '';
    gds.equations = [];
    gds.dim = 0;
    gds.der = [[0 0 0 1 1];zeros(2,5);[1 1 1 0 0]]; 
    gds.jac = '';
    gds.jacp = '';
    gds.hess = '';
    gds.hessp = '';
    gds.tensor3 = '';
    gds.tensor4 = '';
    gds.tensor5 = '';
    gds.point = '';
    gds.type = '';
    gds.discretization.ntst = 20;
    gds.discretization.ncol = 4;
    gds.period = 1;
    gds.plot2 = '';
    gds.plot3 = '';
    gds.PRC='';
    gds.dPRC='';
    gds.open.figuur = 0;
    gds.open.continuer = 0;
    gds.open.numeric_fig = 0;
    gds.open.D2 = 0;
    gds.open.D3 = 0;
    gds.open.PRC = 0;
    gds.open.dPRC = 0;
    gds.open.integrator = 0;
    gds.integrator = [];
    gds.integrator.method = 'ode45';
    gds.integrator.options = [];
    gds.integrator.tspan = [0 1];
    gds.numeric = [];
    gds.numeric.O = {'time' 1;'coordinates' 1;'parameters' 0'};
    gds.numeric.EP = {'coordinates' 1;'parameters' 1;...
        'testfunctions' 0;'eigenvalues' 0;'current stepsize' 0};
    gds.numeric.LC = {'parameters' 1;'period' 1;...
        'testfunctions' 0;'multipliers' 0;'current stepsize' 0;...
        'PRC' 0;'dPRC' 0;'Input' 0};
    gds.numeric.PD = {'parameters' 1;'period' 1;...
        'testfunctions' 0;'multipliers' 0;'current stepsize' 0};
    gds.numeric.Hom = {'parameters' 1;'period' 1;...
        'testfunctions' 0;'multipliers' 0;'current stepsize' 0};
    gds.numeric.HSN = {'parameters' 1;'period' 1;...
        'testfunctions' 0;'multipliers' 0;'current stepsize' 0};
    gds.diagram = 'diagram';
    gds.T = [];
    gds.eps0 = [];
    gds.eps1 = [];
    gds.extravec = [1 0 0];
    gds.t = [];
    gds.epsilon = [];
end

