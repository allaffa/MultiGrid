function [x,y,sol,Rg] = adaptive_mesh_refinement2D(nx0, x0, x1, ny0, y0, y1, epsilon, maxiter, function_string)

    count_iter = 0;

    % n0 = initial value for the number of mesh nodes
    x = linspace(x0, x1, nx0);
    x = x';
    
    % n0 = initial value for the number of mesh nodes
    y = linspace(y0, y1, ny0);
    y = y';    
    
    % Construct and Solve the problem on the coarse mesh
    f = str2func(function_string);
    [A,rhs] = f( x, y );
    sol = A\rhs;
    
    % Compute the approximation of the gradient norm of the error
    Rg = Recovered_gradient(x,y,sol);
    Rg_norm = sum(Rg,2);
    
    display(strcat('Iteration counter: ', num2str(count_iter), ' - Maximal gradient norm: ', num2str(max(Rg_norm))));
    
    % Indentify regiorn of the domain where the gradient norm of the error
    % exceeds the prescribed thhreshold epsilon
    refining_x_indices = find( Rg(:,1) > epsilon );
    refining_y_indices = find( Rg(:,2) > epsilon );
    
    while( size(refining_x_indices,1)>0 && size(refining_y_indices,1)>0  && count_iter<maxiter )
        
        x_indices = [];
        
        for i = 1:size(refining_x_indices,1)
            
            x_indices = [x_indices; mod(refining_x_indices(i),size(x,1))];
            
        end
                
        x_indices = sort(x_indices);
        x_indices = unique(x_indices);
        %x_indices(x_indices==0) = [];
        %x_indices(x_indices==1) = [];
        
        new_x_coordinate = [];

        %Refine the sections with too high gradient norm
        for i = 1:size(x_indices,1)        
                new_x_coordinate = [ new_x_coordinate; ...
                ( x( x_indices(i)-1 ) + x( x_indices(i) ) )/2; ...
                ( x( x_indices(i) ) + x( x_indices(i)+1 ) )/2];
        end
        
        y_indices = [];
       
        for i = 1:size(refining_y_indices,1)
            
            y_indices = [y_indices; mod(refining_y_indices(i),size(y,1))];
            
        end        
       
        y_indices = sort(y_indices);
        y_indices = unique(y_indices); 
        y_indices(y_indices==0) = [];
        y_indices(y_indices==1) = [];
        
        new_y_coordinate = [];

        %Refine the sections with too high gradient norm
        for i = 1:size(y_indices,1)        
                new_y_coordinate = [ new_y_coordinate; ...
                ( y( y_indices(i)-1 ) + y( y_indices(i) ) )/2; ...
                ( y( y_indices(i) ) + y( y_indices(i)+1 ) )/2];
        end        
        
        x = [x; new_x_coordinate];
        x = sort(x);
        x = unique(x);
        
        y = [y; new_y_coordinate];
        y = sort(y);
        y = unique(y);        
        
        count_iter = count_iter + 1;
        
        display('Constructing the problem...');
        [A,rhs] = f( x, y );
        display('Solving the problem...');
        sol = gmres(A,rhs,50,1e-8,200, diag(diag(A))) ;
        %sol = gmres(A,rhs,50,1e-8,1000) ;
    
        % Compute the approximation of the gradient norm of the error
        Rg = Recovered_gradient(x,y,sol);
        Rg_norm = vecnorm(Rg,2,2);

        display(strcat('Iteration counter: ', num2str(count_iter), ' - Maximal gradient norm: ', num2str(max(Rg_norm))));

        % Indentify regiorn of the domain where the gradient norm of the error
        % exceeds the prescribed thhreshold epsilon
        refining_x_indices = find( Rg(:,1) > epsilon );
        refining_y_indices = find( Rg(:,2) > epsilon );

        
    end

return;

end