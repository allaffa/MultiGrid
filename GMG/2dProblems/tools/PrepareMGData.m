function MGData = PrepareMGData( Nx, x1, x2, ratiox, Ny, y1, y2, ratioy, relaxing, problem )    


    MGData = cell(2,1,1);
    nx = (Nx-1)/(ratiox) + 1;
    ny = (Ny-1)/(ratioy) + 1;
    
    N = Nx * Ny;
    n = nx * ny;

    % Fine level
    [A,b] = problem(x1,x2,y1,y2,Nx,Ny);
    MGData{1}.A = A;
    [L,U] = lu(A);
    MGData{1}.L=L;
    MGData{1}.U=U;    
    MGData{1}.P = speye(N);
    MGData{1}.b = b;
    MGData{1}.u = zeros(N,1);
    MGData{1}.R = build_restriction_fullweighting( Nx, nx, x1, x2, Ny, ny, y1, y2 );    
    
    %Coarse level
    [A,~] = problem(x1,x2,y1,y2,nx,ny);
    MGData{2}.A = A;
    [L,U] = lu(A);
    MGData{2}.L=L;
    MGData{2}.U=U;    
    MGData{2}.P = build_interpolant( Nx, nx, x1, x2, Ny, ny, y1, y2 );
    
    MGData{2}.u = zeros(n,1);
    MGData{2}.R = speye(n);
        
    % Construction of smoothers for different levels
    
    if(strcmp(relaxing, 'Jacobi'))
       
        MGData{1}.Relax = sparse( diag(diag(MGData{1}.A)) );
        MGData{2}.Relax = sparse( diag(diag(MGData{2}.A)) );
    
    elseif(strcmp(relaxing, 'GaussSeidel'))    
        
        MGData{1}.Relax = sparse( tril(MGData{1}.A) );
        MGData{2}.Relax = sparse( tril(MGData{2}.A) );
        
    end
    
    return
    
    
end