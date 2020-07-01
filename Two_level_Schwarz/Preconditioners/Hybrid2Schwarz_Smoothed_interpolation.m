function [SP] = Hybrid2Schwarz_Smoothed_interpolation(A, MGData, DD_data, contrib, colouring)
     
    % here colouring is used as a boolean variable

    switch colouring
        
        case 0

            n_domains = size(DD_data,1);

            aux = sparse( size(A,1), size(MGData{2}.P,2) );

            C = A * MGData{2}.P;

            % Parallelizable for loop
            for dom = 1:n_domains       

                 temp = DD_data{dom}.L \ C(DD_data{dom}.global_indices, :);
                 temp = DD_data{dom}.U \ temp;
                 aux(DD_data{dom}.global_indices,:) = aux(DD_data{dom}.global_indices,:) + temp;

            end

            contrib = 1./contrib;
            Diag = spdiags(contrib(:),0,size(A,1),size(A,1));
            aux = Diag * aux;

            SP = MGData{2}.P - aux;
            
            
        case 1
            
 
            n_domains = size(DD_data,1);

            aux = sparse( size(A,1), size(MGData{2}.P,2) );

            C = A * MGData{2}.P;

            % Parallelizable for loop
            for label = 1: 9
                for dom = 1:n_domains       

                    if(DD_data{dom}.colouring == label)
                         temp = DD_data{dom}.L \ C(DD_data{dom}.global_indices, :);
                         temp = DD_data{dom}.U \ temp;
                         aux(DD_data{dom}.global_indices,:) = aux(DD_data{dom}.global_indices,:) + temp;
                    end

                end
            end

            contrib = 1./contrib;
            Diag = spdiags(contrib(:),0,size(A,1),size(A,1));
            aux = Diag * aux;

            SP = MGData{2}.P - aux;
            
           
        
    end
    
end