function [final_mesh, omega] = dynamic_mesh_refinement1D(n0, x0, x1, epsilon, function_string)

    count_iter = 0;

    % n0 = initial value for the number of mesh nodes
    computational_domain = linspace(x0, x1, n0);
    computational_domain = computational_domain';

    physical_domain = computational_domain;    % Construct and Solve the problem on the coarse mesh
    f = str2func(function_string);
    [A,rhs] = f(physical_domain);
    sol = A\rhs;
    
    discrepancy = 1;
    
    while( discrepancy>epsilon )
        
        gradient_solution = Centered_gradient(physical_domain, sol);
        omega = monitor_function1D(physical_domain, gradient_solution); 
        
        temp_physical_domain = zeros(size(computational_domain,1),1);
        temp_physical_domain(1) = x0;
        
        for i=1:size(omega,1)          
            temp_physical_domain(i+1) = temp_physical_domain(i) + 1/omega(i) ;     
        end
        
        temp_physical_domain(end) = temp_physical_domain(end-1) + 1/omega(end-1);
        
        % renormalize 
        temp_physical_domain = x0 + temp_physical_domain / (temp_physical_domain(end) - temp_physical_domain(1)) * (x1-x0);
        
        
        discrepancy = max( physical_domain - temp_physical_domain ); 
        physical_domain = temp_physical_domain; 
        
        [A,rhs] = f(physical_domain);
        sol = A\rhs;
        
        display(strcat('Iteration counter: ', num2str(count_iter), ' - discrepancy: ', num2str(discrepancy)));
         
    end
    
    final_mesh = physical_domain; 

return 
end