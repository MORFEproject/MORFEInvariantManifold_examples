function index = HTHom_nextIndex(settings)


HTHomds = settings.htds;
gds = HTHomds.gds;


number = 0; %number of sparams equal to zero
for i = 1:size(gds.SParams,1)
    if gds.SParams{i,2} == 0
        number = number+1;
    end
end

if number < (size(gds.SParams,1)-1)
    if HTHomds.index == 0 %ConnA_ConnB
        HTHomds.index = 1;

    elseif HTHomds.index == 1 %ConnB_ConnB
    end
    
elseif number == (size(gds.SParams,1)-1)
    if HTHomds.index == 0 %ConnA_ConnC
        HTHomds.index = 2;
    elseif HTHomds.index == 1 %ConnB_ConnC
        HTHomds.index = 2;
    else %HTHomds.index == 2
    end
    
else %ConnC_ConnD of ConnD_ConnD
    

    if HTHomds.index == 2 %ConnC_ConnD
        HTHomds.index = 3;

    else %ConnD_ConnD
    end
    
end
index = HTHomds.index;

end
