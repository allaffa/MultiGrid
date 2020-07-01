function Pr = Hybrid2Schwarz_Preconditioner_modified2(A, SP, MGData, DD_data, contrib, rF, colour)
        
    global DD_time_modified2 Ac_time_modified2 P_time_modified2 Pbar_time_modified2

    n_domains = size(DD_data,1);

    v = zeros(size(A,1),1);
       
    %% This part can be run independently
    
    % Domain Decomposition preconditioner
    % Parallelizable for loop
    
    start = cputime;
    
    switch colour
    
        case 0 
            
            for dom = 1:n_domains       

                z = DD_data{dom}.L \ ( rF(DD_data{dom}.global_indices) ); 

                z = DD_data{dom}.U \ z;

                v(DD_data{dom}.global_indices) = v(DD_data{dom}.global_indices) + z;

            end
            
        case 1
            
            for label = 1: 9 
                
                for dom = 1:n_domains       

                    if (DD_data{dom}.colouring == label)
                        
                        z = DD_data{dom}.L \ ( rF(DD_data{dom}.global_indices) ); 

                        z = DD_data{dom}.U \ z;

                        v(DD_data{dom}.global_indices) = v(DD_data{dom}.global_indices) + z;
                        
                    end

                end
                
            end 
    end
    
    v = contrib.\v;
    
    finish = cputime;
    
    DD_time_modified2 = DD_time_modified2 + (finish - start);
    
    %% This part can be run independently
    w = zeros(size(A,1),1);
    
    % Application of Multigrid preconditioner
    
    start = cputime;   
    wmg = MGData{1}.R * rF;    
    finish = cputime;
    
    P_time_modified2 = P_time_modified2 + (finish - start);    
    
    
    start = cputime;    
    wmg = MGData{2}.L \ wmg;
    wmg = MGData{2}.U \ wmg;   
    finish = cputime;
    
    Ac_time_modified2 = Ac_time_modified2 + (finish-start);
    
    %% Application of smoothed interpolator
    
    start = cputime;
    wmg = SP * wmg;
    finish = cputime;
    
    Pbar_time_modified2 = Pbar_time_modified2 + (finish-start);
    
    %% Gathering operation
    
    Pr = v + wmg;
        
end



