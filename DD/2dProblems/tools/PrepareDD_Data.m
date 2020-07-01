function [DD_data, count_contributions] = PrepareDD_Data( Nx, x1, x2, ndomains_x, Ny, y1, y2, ndomains_y, overlap, problem )    

    % The labeling of the subdomains is as follows


    %   (ndomains_y - 1)*ndomains_x + 1                                                        ...  ndomains_y *  ndomains_x
    %
    %   .
    %   .
    %   .
    %
    %   2 * ndomains_x + 1
    %
    %   ndomains_x     + 1
    %
    %   1                                           2               3               ...             ndomains_x



    [Ag,rhsg] = problem ( x1,x2,y1,y2,Nx,Ny );

    % Total number of subdomains
    n_domains = ndomains_x * ndomains_y;

    % size of each subdomain in the x-direction
    sizex_d = floor( Nx / ndomains_x );

    % size of each subdomain in the y-direction
    sizey_d = floor( Ny / ndomains_y );
    
    %total size of a subdomain
    size_d = sizex_d * sizey_d;
    
    % Data structure to save quantities of each subdomain 
    DD_data = cell(n_domains, 1);
    
    % vector to employ the number of contributions added in the node
    count_contributions = zeros(Nx*Ny,1);


    % Initial partition into subdomains without overlap

    for domy = 1:ndomains_y
    
        for domx = 1:ndomains_x

            % computaiton of the subdomain index
            n = (domy-1)*ndomains_x + domx;
            
            lnodes = zeros(Nx*Ny,1); 
            
            if (domy ~= ndomains_y) && (domx ~= ndomains_x)
                
                for k = 1:sizey_d
                    start = (domy-1) * sizey_d * Nx + (domx-1) * sizex_d + 1 + (k-1) * Nx;
                    finish = (domy-1) * sizey_d * Nx + domx * sizex_d + (k-1) * Nx;
                    lnodes( start : finish ) = ones( finish - start + 1 , 1 );
                end
                
            elseif (domy == ndomains_y) && (domx ~= ndomains_x)   
                
                for k = (domy-1)*sizey_d + 1: Ny
                    start = (k-1) * Nx + (domx-1) * sizex_d + 1;
                    finish = (k-1) * Nx + domx * sizex_d;
                    lnodes( start : finish ) = ones( finish - start + 1 , 1 );
                end
                
            elseif (domy ~= ndomains_y) && (domx == ndomains_x)   
                
                for k = 1:sizey_d
                    start = (domy-1) * sizey_d * Nx + (domx-1) * sizex_d + 1 + (k-1) * Nx;
                    finish = (domy-1) * sizey_d * Nx + k * Nx;
                    lnodes( start : finish ) = ones( finish - start + 1 , 1 );
                end
                
            else 
                
               for k = (domy-1)*sizey_d + 1: Ny
                    start = (k-1) * Nx + (domx-1) * sizex_d + 1;
                    finish = k * Nx;
                    lnodes( start : finish ) = ones( finish - start + 1 , 1 );
                end                
                
            end

                
            lnodes_old = lnodes;       

            if (domy == 1 && domx == 1) || (domy == 1 && domx == ndomains_x) || (domy == ndomains_y && domx == 1) || (domy == ndomains_y && domx == ndomains_x)

                while( nnz(lnodes) < (1+2*overlap)*size_d )

                    lnodes = abs(Ag) * lnodes;
                    lnodes = spones( lnodes );

                end

            elseif  (domy == 1 && domx ~= 1) || (domy == 1 && domx ~= ndomains_x) || (domy == ndomains_y && domx ~= 1) || (domy == ndomains_y && domx ~= ndomains_x)

               while( nnz(lnodes) < (1+3*overlap)*size_d )

                    lnodes = abs(Ag) * lnodes;
                    lnodes = spones( lnodes );

                end            

            else 

               while( nnz(lnodes) < (1+4*overlap)*size_d )

                    lnodes = abs(Ag) * lnodes;
                    lnodes = spones( lnodes );

               end

            end
            
            count_contributions = count_contributions + spones( lnodes );
            
            R = sparse(nnz(lnodes), Nx*Ny);
        
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
            
            % Keep track of the overlap

            DD_data{n}.R              = R;
            DD_data{n}.A              = A;
            DD_data{n}.L              = L;
            DD_data{n}.U              = U;            
            DD_data{n}.rhs            = rhs;
            DD_data{n}.global_indices = indices; 
            DD_data{n}.overlap        = find(lnodes - lnodes_old);
            
        end

    end

    
    
end