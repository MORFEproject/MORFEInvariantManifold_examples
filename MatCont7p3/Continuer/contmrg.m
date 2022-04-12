function res = contmrg(options, opt)
%
% res = CONTMRG(options, opt)
%
% Merges options in 'opt' to 'options'
% Options in 'options' will be overwritten by defined options in 'opt'

if isempty(opt)
  res = options;
  return;
end

allopt = contidx;

[m,n] = size(allopt);

for i=1:m
    if isfield(opt, strtrim(allopt(i, :)))
        
        eval(['vals = opt.' allopt(i,:) ';']);
    else
        eval(['vals = options.' allopt(i,:) ';']);
        
    end
    
    if strmatch('IgnoreSingularity',allopt(i,:))
        eval(['vals1 = options.' allopt(i,:) ';']);
        values = mergearray(vals1,vals);
        eval(['options.' allopt(i,:) '= values;']);
    elseif ~isempty(vals)
        eval(['options.' allopt(i,:) '= vals;']);
    end
    
end

res = options;

function vals1 = mergearray(vals1,val)
    for i=1:length(val)
        if isempty(vals1)
            str='[]';
            st = [];
        else
            str = mat2str(vals1);
            st  = vals1;
        end         
        x = find(st==val(i)) ;
        if ~isempty(x)
            continue
        end
        if (size(vals1,2)==1)
            vals1 = str2num(str2mat(sprintf('[%s%c%s]',str,char(32),num2str(val(i)))));
        else      
            str(end)='';
            vals1 = str2num(str2mat(sprintf('%s%c%s]',str,char(32),num2str(val(i)))));
            
        end
    end
    if all(size(vals1)==[1 0])
        vals1=[];
    end
