function index = HTHet_nextIndex(settings)

HTHetds = settings.htds;
gds = HTHetds.gds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START ORIGINAL CODE %%%%%%%%%%%%%%%%%%%%%%%%


number = 0; %number of sparams equal to zero
for i = 1:size(gds.SParams,1)
    if gds.SParams{i,2} == 0
        number = number+1;
    end
end

dim_npos = size(gds.UParams,1);

   
if number < (dim_npos-1)
    if HTHetds.index == 0 %ConnA_ConnB      
        HTHetds.index = 1;
    elseif HTHetds.index == 1 %ConnB_ConnB                     
    end
    
elseif (size(gds.SParams,1)-number) > 0        
           
    if HTHetds.index == 0 %ConnA_ConnC            
        HTHetds.index = 2;
    elseif HTHetds.index == 1 %ConnB_ConnC     
        HTHetds.index = 2;
    else %HTHetds.index == 2       
    end           
    
else %ConnC_ConnD of ConnD_ConnD
    if HTHetds.index == 1
        HTHetds.index = 3;
    elseif HTHetds.index == 2 %ConnC_ConnD
        HTHetds.index = 3;
    else %ConnD_ConnD 
    end
    
end                       

index = HTHetds.index;
end
