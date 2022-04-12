function [equations, co, pa] = rsym(equations, co, pa)
    if ~isempty(co)
        co = strsplit(co, ',');
    else
        co = {};
    end
    if ~isempty(pa)
        pa = strsplit(pa, ',');
    else
        pa = {};
    end

    renamings_co = [co; cellfun(@(s) strcat('cor_', s), co,'UniformOutput',0)];
    renamings_pa = [pa; cellfun(@(s) strcat('par_', s), pa,'UniformOutput',0)];
    
    filterspaces = @(z) z(z ~= ' ');
    
    eqn = cell(1, size(equations, 1));
    for row = 1:size(equations, 1)
        eqn{row} = filterspaces(equations(row, :));
        %maxlen = max(maxlen, length(eqn{row}));
    end
    for rename = [renamings_co, renamings_pa]
        eqn = regexprep(eqn, strcat('(^|\W)', rename{1}, '(\W|$)'), strcat('$1', rename{2}, '$2'));
        eqn = regexprep(eqn, strcat('(^|\W)', rename{1}, '(\W|$)'), strcat('$1', rename{2}, '$2')); %pattern overlap
    end
    maxlen = 0;
    for row = 1:size(equations, 1)
        maxlen = max(maxlen, length(eqn{row}));
    end    
    equations = repmat(' ', size(equations, 1), maxlen);
    for row = 1:size(equations, 1)
        equations(row, :) = [eqn{row},   repmat(' ', 1, maxlen - length(eqn{row}))];
    end
    if ~isempty(renamings_co)
        co = strjoin(renamings_co(2,:), ',');
    end
    if ~isempty(renamings_pa)
        pa = strjoin(renamings_pa(2,:), ',');
    end
    
end

