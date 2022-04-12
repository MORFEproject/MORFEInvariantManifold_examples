function gds_struct = gds_from_sbml(model)
%GDS_FROM_SBML Using a SBML_model structure, return MATCONT gds struct

if (~isSBML_Model(model))
    error('create_from_sbml(): argument must be an SBML Model');
end

if (~isempty(model.id))
    model_name = model.id;
elseif (~isempty(model.name))
    model_name = model.name;
else
    error('gds_from_sbml(): SBML Model must have a name');
end

parameters_names = [];
parameters_values = [];
assignments_names = [];
assignments_formulas = [];
% collect all the parameters from the reactions
for reaction = model.reaction
    [names values] = GetParameterFromReaction(reaction);
    for i = 1:length(names)
        parameters_names = [parameters_names names(i)];
        parameters_values = [parameters_values values(i)];
    end
end
% collect the global parameters
for i=1:length(model.parameter)
    parameter = model.parameter(i);
    if parameter.constant
        parameters_names = [parameters_names cellstr(parameter.id)];
        parameters_values = [parameters_values parameter.value];
    else
        assignments_names = [assignments_names cellstr(parameter.id)];
        for rule = model.rule
            if strcmp(rule.variable,parameter.id)
                formula = substitute_functions(rule.formula,model);
            end
        end
        assignments_formulas = [assignments_formulas cellstr(formula)];
    end
end
waitbar(0.1);
species = AnalyseSpecies(model);
waitbar(0.2);
n_species = length(species);
reactions = model.reaction;
compartments = model.compartment;
if (~isempty(model.time_symbol))
    time_symbol = model.time_symbol;
else
    time_symbol = 't';
end
waitbar(0.3);

% MATCONT uses coordinates, SBML uses species
coordinates_names = [];
coordinates_string = '[';
coordinates_equations = [];
coordinates_initial_values = [];
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
    waitbar(0.3+0.2*i/n_species);
end
n_assignments = length(assignments_names);
n_parameters = length(parameters_names);
coordinates_string = strcat(coordinates_string,']');
n_coordinates = length(coordinates_names);
assert(n_coordinates > 0, 'gds_from_sbml(): Model has to have coordinates');
% ignore the compartments, they have a 'size' set, important for matcont???
for compartment = compartments
    compartment_name = compartment.id;
    coordinates_equations = strrep(coordinates_equations,...
        compartment_name,'1.0');
end
waitbar(0.7);
% sbml adds the name of the reaction to the
% parameter if there might be overlap, but in
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
waitbar(0.8);

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
    syms(char(parameters_names(i)));
end
gds_struct.time = {char(time_symbol),0};
for i = 1:n_assignments
    lhs=char(assignments_names(i));
    rhs=char(assignments_formulas(i));
    rhs=substitute_functions(rhs,model);
    rhs=replace_nthroot_by_power(rhs);
    gds_struct.equations = strvcat(gds_struct.equations,strcat(lhs,'=',rhs));
    syms(lhs);
end
for i = 1:n_coordinates
    lhs = [char(coordinates_names(i)) ''''];
    equation = char(coordinates_equations(i));
    equation = substitute_functions(equation,model);
    rhs = char(eval(char(equation)));
    %rhs = substitute_assignments(rhs,...
    %    assignments_names,assignments_formulas);
    rhs = substitute_functions(rhs,model);
    rhs = replace_nthroot_by_power(rhs);
    %rhs = substitute_assignments(rhs,...
    %    assignments_names,assignments_formulas);
    lhs = strtrim(lhs);
    rhs = strtrim(rhs);
    gds_struct.equations = strvcat(gds_struct.equations,[lhs '=' rhs]);
end
gds_struct.dim = n_coordinates;
waitbar(1);

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

function formula = substitute_functions(formula, model)
subs_applied = true;
prev_formula = formula;
while subs_applied
    subs_applied = false;
    for function_definition = model.functionDefinition
        formula = SubstituteFunction(formula,function_definition);
    end
    subs_applied = ~strcmp(formula,prev_formula);
    prev_formula = formula;
end
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
        exp = strcat('\<',assignments_names(i),'\>');
        repstr = assignments_formulas(i);
        formula = regexprep(str,exp,repstr);
        rule_applied = rule_applied + strcmp(str, formula)==false;
    end
    iterations_left = iterations_left - 1;
end
assert(rule_applied == 0, ...
    'substitute_assignments(): Cyclic dependency of rules dedected');
end

function formula = replace_nthroot_by_power(formula)
match='\<nthroot\s*(';
replacement='power(';
pos = regexp(formula,match,'once');
while ~isempty(pos)
    formula = regexprep(formula,match,replacement,'once');
    formula = [formula ' '];% allow to go one to far
    i=pos+length(replacement);
    bracket_count = 1;
    ch = formula(i);
    while bracket_count > 0
        switch ch
            case '('
                bracket_count = bracket_count + 1;
            case ')'
                bracket_count = bracket_count - 1;
            case ','
                if bracket_count == 1
                    begin_second_arg = i+1;
                end
        end
        i = i+1;
        ch = formula(i);
    end
    end_second_arg = i-1;
    second_argument = formula(begin_second_arg:end_second_arg);
    formula = strcat(formula(1:begin_second_arg-1),...
        '1/',second_argument,formula(end_second_arg+1:end));
    pos = regexp(formula,match,'once');
end
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
    gds.der = [[1 1 1 1 1];zeros(3,5);]; 
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

