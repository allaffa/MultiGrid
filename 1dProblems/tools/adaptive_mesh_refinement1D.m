function [final_mesh, Rg] = adaptive_mesh_refinement1D(n0, x0, x1, epsilon, function_string)

    count_iter = 0;

    % n0 = initial value for the number of mesh nodes
    x = linspace(x0, x1, n0);
    x = x';
    
    % Construct and Solve the problem on the coarse mesh
    f = str2func(function_string);
    [A,rhs] = f(x);
    sol = A\rhs;
    
    % Compute the approximation of the gradient norm of the error
    Rg = Recovered_gradient(x,sol);
    
    display(strcat('Iteration counter: ', num2str(count_iter), ' - Maximal gradient norm: ', num2str(max(Rg))));
    
    % Indentify regiorn of the domain where the gradient norm of the error
    % exceeds the prescribed thhreshold epsilon
    % In need to add an offest of 1 because Rg is only defined at internal
    % points
    refining_indices = find( Rg > epsilon ) + 1;
    
    while( size(refining_indices,1)>0 )
        
        new_node_coordinate = [];
        
        %Refine the sections with too high gradient norm
        for i = 1:size(refining_indices,1)
            
            new_node_coordinate = [ new_node_coordinate; ...
                ( x( refining_indices(i)-1 ) + x( refining_indices(i) ) )/2; ...
                ( x( refining_indices(i) ) + x( refining_indices(i)+1 ) )/2];
            
        end
        
        x = [x; new_node_coordinate];
        x = sort(x);
        x = unique(x);
        
        count_iter = count_iter + 1;
        
        [A,rhs] = fd1d_bvp_test03 (x);
        sol = A\rhs;
    
        % Compute the approximatiosn of the gradient norm of the error
        Rg = Recovered_gradient(x,sol);
        refining_indices = find( Rg > epsilon ) + 1;
        
        display(strcat('Iteration counter: ', num2str(count_iter), ' - Maximal gradient norm: ', num2str(max(Rg))));
        
    end

    final_mesh = x;

return;

end