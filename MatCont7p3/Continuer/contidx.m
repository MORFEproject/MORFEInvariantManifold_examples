function opt = contidx(option)
%
% opt = CONTIDX(option)
%
% Returns index of an option to be used in option vector
%
% Input : option: string with option name
% Output: index : index of option in option vector
%         maxopt: number of supported options

Names = [
    'InitStepsize     '
    'MinStepsize      '
    'MaxStepsize      '
    'MaxCorrIters     '
    'MaxNewtonIters   '
    'MaxTestIters     '
    'MoorePenrose     '
    'SymDerivative    '
    'SymDerivativeP   '
    'Increment        '
    'FunTolerance     '
    'VarTolerance     '
    'TestTolerance    '
    'Singularities    '
    'MaxNumPoints     '
    'Backward         '
    'CheckClosed      '
    'TestFunctions    '
    'WorkSpace        '
    'Locators         '
    'Adapt            '
    'IgnoreSingularity'
    'ActiveParams     '
    'Multipliers      '
    'Eigenvalues      '
    'Userfunctions    '
    'UserfunctionsInfo'
    'PRC              '         
    'dPRC             '         
    'Input            '    
    'ActiveUParams    '
    'ActiveSParams    '
    'ActiveSParam     ' 
    'TSearchOrder     '
%     'ActiveT          '
    ];

   
if nargin==0
  opt = Names;
  return;
end

[m,n] = size(Names);
names = lower(Names);

option = deblank(lower(option));

index = 0;
for i=1:m
  if strcmp(option, deblank(names(i,:)))
    index = i;
  end
end

if index ~= 0
  opt = Names(index,:);
else
  opt = [];
end

%SD:option to index nr

