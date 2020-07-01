function Pr = Hybrid2Schwarz_Preconditioner_modified(A, MGData, DD_data, contrib, rF)
 
global DD_time_modified Ac_time_modified P_time_modified 

    n_domains = size(DD_data,1);
    
    v = zeros(size(A,1),1);
    
    %% This part can be run independently
    
    % Domain Decomposition preconditioner
    % Parallelizable for loop
    start = cputime;
    for dom = 1:n_domains       
        
        z = DD_data{dom}.L \ ( rF(DD_data{dom}.global_indices) ); 
        
        z = DD_data{dom}.U \ z;
        
        v(DD_data{dom}.global_indices) = v(DD_data{dom}.global_indices) + z;
        
    end
    
    v = contrib.\v;
    
    finish = cputime;
    
    DD_time_modified = DD_time_modified + (finish - start);
    
    %% This part can be run independently
    w = zeros(size(A,1),1);
    
    % Application of Multigrid preconditioner
    start = cputime;
    wmg = MGData{1}.R * rF;
    finish = cputime;
    
    P_time_modified = P_time_modified + (finish - start);
    
    start = cputime;
    wmg = MGData{2}.L \ wmg;
    wmg = MGData{2}.U \ wmg;
    finish = cputime;
    
    Ac_time_modified = Ac_time_modified + (finish - start);
    
    start = cputime;
    wmg = MGData{2}.P * wmg;
    finish = cputime;
    
    P_time_modified = P_time_modified + (finish - start);
    
    % Premultiplication by Af before applying Domain Decomposition
    % preconditioner
    
    start = cputime;
    Awmg = A * wmg;
    
    % Domain Decomposition preconditioner
    % Parallelizable for loop
    for dom = 1 : n_domains
        
      wp = DD_data{dom}.L \  Awmg(DD_data{dom}.global_indices);    
      
      wp = DD_data{dom}.U \ wp;
      
      w(DD_data{dom}.global_indices) = w(DD_data{dom}.global_indices) + wp;
    
    end
    
    w = wmg - contrib.\w;
    finish = cputime;
    
    DD_time_modified = DD_time_modified + (finish-start);
    
    %% Gathering operation
    
    Pr = v + w;
        
end



