function value = Substitute(original_formula, model)
%SUBSTITUE see documentation original Substitute in 
%     SBMLToolbox_*_src/toolbox/Convenience/Substitute.m

formula = original_formula;

% handle easy case
value = str2double(original_formula);
if ~isnan(value)
    return;
end

% put everything in MATLAB and evaluate the formule
% (I'm not sure if all cases are handled, for instance if size is not set
% of compartment, should then be the volume taken?)
if (model.SBML_level > 1 && ~isempty(model.time_symbol))
    assert(exist(model.time_symbol,'var')==false);
    assignin('caller',model.time_symbol,0);
end
for specie = model.species
    if specie.isSetInitialAmount
        assert(exist(specie.id,'var')==false);
        assignin('caller',specie.id,specie.initialAmount);
        %exp = strcat('\<',specie.id,'\>');
        %repstr = specie.name;
        %formula = regexprep(formula,exp,repstr);
    end
end
for parameter = model.parameter
    if parameter.isSetValue
        assert(exist(parameter.id,'var')==false);
        assignin('caller',parameter.id,parameter.value);
        %exp = strcat('\<',parameter.id,'\>');
        %repstr = parameter.name;
        %formula = regexprep(formula,exp,repstr);
    end
end
for compartment = model.compartment
    if compartment.isSetSize
        assert(exist(compartment.id,'var')==false);
        assignin('caller',compartment.id,compartment.size);
        %exp = strcat('\<',compartment.id,'\>');
        %repstr = compartment.name;
        %formula = regexprep(formula,exp,repstr);
    end
end

% this replaces all rules in the original formula
rule_applied = 1;
iterations_left = length(model.rule) + 1;
while rule_applied > 0 && iterations_left > 0
    rule_applied = 0;
    for rule = model.rule
        str = formula;
        exp = strcat('\<',rule.variable,'\>');
        repstr = rule.formula;
        formula = regexprep(str,exp,repstr);
        rule_applied = rule_applied + strcmp(str, formula)==false;
    end
    iterations_left = iterations_left - 1;
end
assert(rule_applied == 0, ...
    'Substitute(): Cyclic dependency of rules dedected');

% this substitutes all functions in the original formula
subs_applied = true;
iterations_left = length(model.functionDefinition);
prev_formula = formula;
while subs_applied && iterations_left > 0
    subs_applied = false;
    for function_definition = model.functionDefinition
        formula = SubstituteFunction(formula,function_definition);
    end
    subs_applied = ~strcmp(formula,prev_formula);
    prev_formula = formula;
    iterations_left = iterations_left - 1;
end
assert(subs_applied == 0, ...
    'Substitute(): Cyclic dependency of function definitions dedected');

try
    value = evalin('caller',formula);
catch
    error('Substitute(): Ill formed formula');
end

end

