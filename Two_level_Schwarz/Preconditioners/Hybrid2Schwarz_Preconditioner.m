function Pr = Hybrid2Schwarz_Preconditioner(A, MGData, DD_data, contrib, rF)

global DD_time Ac_time P_time 

    %Restriciton of the vector to the coarse level
    start = cputime;
    rC = MGData{1}.R * rF;
    finish = cputime;
    
    P_time = P_time + (finish - start);

    % Coarse level preconditioning
    start = cputime;
    PrC = MGData{2}.L \ rC;
    PrC = MGData{2}.U \ PrC;
    finish = cputime;
    
    Ac_time = Ac_time + (finish - start);

    %Prolungation of the coarse level preconditioning
    start = cputime;
    PrF = MGData{2}.P * PrC;
    finish = cputime;
    
    P_time = P_time + (finish - start);
    
    v = zeros(size(A,1),1);
    
    n_domains = size(DD_data,1);
        
    start = cputime;
    aux = ( rF - A * PrF); 
   
    % Parallelizable for loop
    for dom = 1 : n_domains
        w = DD_data{dom}.L \ aux( DD_data{dom}.global_indices );
        w = DD_data{dom}.U \ w;
        v( DD_data{dom}.global_indices ) = v( DD_data{dom}.global_indices ) +  w;
    end
    
    v = contrib.\v;
    
    finish = cputime;
    
    DD_time = DD_time + (finish-start);
    
    % Gathering operation
    Pr = PrF + v;
    
        
end



