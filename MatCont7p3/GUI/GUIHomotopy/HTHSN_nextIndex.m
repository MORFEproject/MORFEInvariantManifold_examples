function index = HTHSN_nextIndex(settings)



HTHSNds = settings.htds;
gds = HTHSNds.gds;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START ORIGINAL CODE %%%%%%%%%%%%%%%%%%%%%%%%



number = 0; %number of sparams equal to zero
for i = 1:size(gds.SParams,1)
    if gds.SParams{i,2} == 0
        number = number+1;
    end
end

if (size(gds.SParams,1)-number) > 0
    if HTHSNds.index == 0
        HTHSNds.index = 1;
    elseif HTHSNds.index == 1
    end
else
    
    if HTHSNds.index == 0
        HTHSNds.index = 2;
    elseif HTHSNds.index == 1
        HTHSNds.index = 2;
    else
    end
    
end

index = HTHSNds.index;
end