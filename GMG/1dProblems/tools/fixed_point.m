function u = fixed_point( N, x1,x2,nu1, relaxing )

    num_levels = 2;
    num_cycles = 1;

    MGData = cell(num_levels,1,1);
    n = (N-1)/(2) + 1;

    % Fine level
    [A,b] = fd1d_bvp_test01 ( N,x1,x2 );
    MGData{1}.A = A;
    MGData{1}.P = speye(N);
    MGData{1}.b = b;
    MGData{1}.u = zeros(N,1);
    MGData{1}.R = build_restriction_injection( N, n, x1, x2 );    
    
        
    % Construction of smoothers for different levels
    
    if(strcmp(relaxing, 'Jacobi'))
        omega = 1;
        MGData{1}.Relax = sparse( diag(diag(MGData{1}.A)) );
    
    elseif(strcmp(relaxing, 'GaussSeidel'))    
        
        
    end
    
    for l = 1:nu1
        
            z = MGData{1}.A * MGData{1}.u;
            MGData{1}.u = (1-omega) * MGData{1}.u  +  ...
                omega * ( MGData{1}.u - MGData{1}.Relax \ z + MGData{1}.Relax \ MGData{1}.b );
        
    end
        
    u = MGData{1}.u;
    
return 
end