function [DD_data, count_contributions] = PrepareDD_Data( N, x1,x2, n_domains, overlap, problem )    

    [Ag,rhsg] = problem ( N,x1,x2 );

    % size of each subdomain 
    size_d = floor( N / n_domains );
    
    % Data structure to save quantities of each subdomain 
    DD_data = cell(n_domains, 1);
    
     % vector to employ the number of contributions added in the node
    count_contributions = zeros(N,1);   
    
    for n = 1 : n_domains-1
        
        % detection of subdomain including overlapping nodes
        
        lnodes = zeros(N,1);      
        lnodes( ((n - 1)*size_d + 1) : (n - 1) *size_d + size_d ) = ones( size_d, 1 );
        
        lnodes_old = lnodes;
        
        if n == 1 
            
            while( nnz(lnodes) < (1+overlap)*size_d )

                lnodes = abs(Ag) * lnodes;
                lnodes = spones( lnodes );

            end
            
        else 
            
           while( nnz(lnodes) < (1+2*overlap)*size_d )

                lnodes = abs(Ag) * lnodes;
                lnodes = spones( lnodes );

            end            
            
        end
        
        
        R = sparse(nnz(lnodes), N);
        
        indices = find(lnodes);
        
        count_contributions = count_contributions + spones( lnodes );
        
        for i = 1 : nnz(lnodes)
           
            R( i, indices(i)) = 1;
            
        end
        
        % Construction of the local stiffness matrix
        A = R * Ag * R';

        % Construction of LU factors for the local stiffness matrix
        [L,U] = lu(A);

        
        % Construction of the local right hand side
        rhs = R * rhsg;
        
        % Keep track of the overlap
        
        DD_data{n}.R              = R;
        DD_data{n}.A              = A;
        DD_data{n}.L              = L;
        DD_data{n}.U              = U;       
        DD_data{n}.rhs            = rhs;
        DD_data{n}.global_indices = indices; 
        DD_data{n}.overlap        = find(lnodes - lnodes_old);
        
        
    end
    
    % the last domain inglobes extra nodes to fit the entire domain
    lnodes = zeros(N,1);      
    lnodes( ((n_domains - 1)*size_d + 1) : N ) = ones( N - ((n_domains - 1)*size_d) , 1 );
    
    lnodes_old = lnodes;

    while( nnz(lnodes) < (1+overlap)*size_d )

        lnodes = abs(Ag) * lnodes;
        lnodes = spones( lnodes );

    end  
    
    R = sparse(nnz(lnodes), N);
    
    count_contributions = count_contributions + spones( lnodes );

    indices = find(lnodes);

    for i = 1 : nnz(lnodes)

        R( i, indices(i) ) = 1;

    end

    % Construction of the local stiffness matrix
    A = R * Ag * R';
    
    % Construction of LU factors for the local stiffness matrix
    [L,U] = lu(A);

    % Construction of the local right hand side
    rhs = R * rhsg;

    DD_data{n_domains}.R              = R;
    DD_data{n_domains}.A              = A;
    DD_data{n_domains}.L              = L;
    DD_data{n_domains}.U              = U;
    DD_data{n_domains}.rhs            = rhs;  
    DD_data{n_domains}.global_indices = indices; 
    DD_data{n_domains}.overlap        = find(lnodes - lnodes_old);
    
    
end