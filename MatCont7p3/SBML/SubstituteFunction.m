function Formula = SubstituteFunction(OriginalFormula, SBMLFunctionDefinition)
% SubstituteFunction 
%       takes 
%           1) a string representation of a formula 
%           2) the SBMLFunctionDefinition structure defining the formula
%       and returns 
%           the formula with the function substituted
%
%
%   EXAMPLE:
%           fD = SBMLFunmctionDefinition 
%               with id = 'g' and math = 'lambda(x,x+0.5)' 
%
%           formula = SubstituteFormula('g(y)', fD)
%           
%                   = 'y+0.5'
%
%           formula = SubstituteFormula('h(y)', fD)
%           
%                   = ''


%  Filename    :   SubstituteFunction.m
%  Description : 
%  Author(s)   :   SBML Development Group <sbml-team@caltech.edu>
%  Organization:   University of Hertfordshire STRI
%  Created     :   11-Feb-2005
%  Revision    :   $Id: SubstituteFunction.m 10822 2010-01-25 16:40:37Z sarahkeating $
%  Source      :   $Source v $
%
%<!---------------------------------------------------------------------------
% This file is part of SBMLToolbox.  Please visit http://sbml.org for more
% information about SBML, and the latest version of SBMLToolbox.
%
% Copyright 2005-2007 California Institute of Technology.
% Copyright 2002-2005 California Institute of Technology and
%                     Japan Science and Technology Corporation.
% 
% This library is free software; you can redistribute it and/or modify it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation.  A copy of the license agreement is provided
% in the file named "LICENSE.txt" included with this software distribution.
% and also available online as http://sbml.org/software/sbmltoolbox/license.html
%----------------------------------------------------------------------- -->


Formula = OriginalFormula;


%check arguments are appropriate
if (~isstruct(SBMLFunctionDefinition))
  error(sprintf('%s', ...
    'first argument must be an SBML functionDefinition structure'));
end;
 
[sbmlLevel, sbmlVersion] = GetLevelVersion(SBMLFunctionDefinition);

if (~ischar(Formula))
    error('SubstituteFunction(OriginalFormula, SBMLFunctionDefinition)\n%s', 'first argument must be a character array containing the id of the function definition');
elseif (~isSBML_FunctionDefinition(SBMLFunctionDefinition, sbmlLevel, sbmlVersion))
    error('SubstituteFunction(OriginalFormula, SBMLFunctionDefinition)\n%s', 'second argument must be an SBML function definition structure');
end;

Formula = LoseWhiteSpace(Formula);

% strfind will find all ocurencies of the function definition
startPoint = strfind(Formula, SBMLFunctionDefinition.id);
if (isempty(startPoint))
  %Formula = '';
  return;
else
  % for each function in the formula, check if use is correct: bracket
  % follows the function name
  for j = startPoint
    funcName = '';
    c = Formula(j);
    while (~strcmp(c, '('))
      funcName = strcat(funcName, c);
      j = j + 1;
      c = Formula(j);
    end;
    
    if (~isequal(funcName, SBMLFunctionDefinition.id))
      Formula = '';
      return;   % FIXME: better give user error
    end;
  end
end;

ElementsOfFuncDef = GetArgumentsFromLambdaFunction(SBMLFunctionDefinition.math);

% get the arguments of the application of the formula
Starts = findstr(Formula, SBMLFunctionDefinition.id);
while ~isempty(Starts)
  StartFunctionInFormula = Starts(1);

  % Get the start of the argument list in the formula, read away one
  % bracket
  j = StartFunctionInFormula + length(SBMLFunctionDefinition.id) + 1;
  c = Formula(j);
  element = '';
  NoElements = 1;
  ElementsInFormula = {};
  % one bracket read, so initialize
  bracket_count = 1;
  % go through the argument list, the argument list stops on it's last
  % bracket
  while (bracket_count>0)
      c = Formula(j);
      if (strcmp(c, ',') && bracket_count==1)
          ElementsInFormula{NoElements} = element;
          element = '';
          NoElements = NoElements + 1;
      elseif (strcmp(c,'('))
          bracket_count = bracket_count+1;
          element = strcat(element,c);
      elseif strcmp(c,')')
          bracket_count = bracket_count-1;
          element = strcat(element,c);
      else
          element = strcat(element, c);
      end;
      j = j + 1;
  end;
  ElementsInFormula{NoElements} = element(1:end-1);
  OriginalFunction = '';
  for i = StartFunctionInFormula:j-1
      OriginalFunction = strcat(OriginalFunction, Formula(i));
  end;

  % check got right number
  if (NoElements ~= length(ElementsOfFuncDef) - 1)
      error('SubstituteFunction(OriginalFormula, SBMLFunctionDefinition)\n%s', 'mismatch in number of arguments between formula and function');
  end;

  % check that same arguments have not been used

  for i = 1:NoElements
    for j = 1:NoElements
      if (strcmp(ElementsInFormula{i}, ElementsOfFuncDef{j}))
        newElem = strcat(ElementsInFormula{i}, '_new');
        ElementsOfFuncDef{j} = newElem;
        ElementsOfFuncDef{end} = strrep(ElementsOfFuncDef{end}, ElementsInFormula{i}, newElem);
      end;
    end;
  end;
  % replace the arguments in function definition with those in the formula
  FuncFormula = '(';
  FuncFormula = strcat(FuncFormula, ElementsOfFuncDef{end});
  FuncFormula = strcat(FuncFormula, ')');
  for i = 1:NoElements
      FuncFormula = strrep(FuncFormula, ElementsOfFuncDef{i}, ElementsInFormula{i});
  end;

  Formula = strrep(Formula, OriginalFunction, FuncFormula);
  Starts = findstr(Formula, SBMLFunctionDefinition.id);
end