classdef InputRestrictions

    properties(Constant)
        INT_g1  = InputRestriction(@(x) isscalar(x) && (floor(x) == x) && (x > 1) , @(x) floor(x), @(x) num2str(x , '%.15g'), 'integer must be > 1');
        INT_g0  = InputRestriction(@(x) isscalar(x) && (floor(x) == x) && (x > 0) , @(x) floor(x), @(x) num2str(x , '%.15g'), 'integer must be > 0');
        INT_ge0 = InputRestriction(@(x) isscalar(x) && (floor(x) == x) && (x >= 0), @(x) floor(x), @(x) num2str(x , '%.15g'), 'integer must be >= 0');
        INT     = InputRestriction(@(x) isscalar(x) && (floor(x) == x)            , @(x) x       , @(x) num2str(x , '%.15g'), 'value is not an integer');
        POS     = InputRestriction(@(x) isscalar(x) && (x > 0)                    , @(x) x       , @(x) num2str(x , '%.15g'), 'number must be > 0');
        POS_0   = InputRestriction(@(x) isscalar(x) && (x >= 0)                   , @(x) x       , @(x) num2str(x , '%.15g'), 'number must be >= 0');
        NUM     = InputRestriction(@(x) isscalar(x)                               , @(x) x       , @(x) num2str(x , '%.15g'), 'must be a number');
        VECTOR  = InputRestriction(@(x) isvector(x) || isempty(x)                 , @(x) x       , @(x) InputRestriction.vector2string(x), 'input must be a vector');
        BOOL    = InputRestriction(@(x) isscalar(x) && ((x == 0) | (x == 1))      , @(x) x~=0    , @(x) InputRestriction.bool2string(x), 'input must be boolean: 0/1 or false/true', @(p, s, sn, varargin) GUISettingSwitch('checkbox', p, s, sn, varargin), 'BOOL' );
        NONE     = InputRestriction(@(x) 1, @(x) x, @(x) x, '');
    end
        
    methods(Static)
        function IR = vector(n)
            IR = InputRestriction(@(x) isvector(x) && length(x) == n , @(x) x       , @(x) InputRestriction.vector2string(x), ['input must be a vector with length ' num2str(n)]);
        end
        
        function IR = vectorOrScalar(n)
            IR = InputRestriction(@(x) (isvector(x) && length(x) == n) || isscalar(x) , @(x) x , @(x) vectorscal2string(x), ['input must be a scalar or a vector with length ' num2str(n)]);
        end        
        
        
    end

end


function s = vectorscal2string(x)
    if isscalar(x)
        s = num2str(x , '%.15g');
    else
        s = InputRestriction.vector2string(x);
    end
end
