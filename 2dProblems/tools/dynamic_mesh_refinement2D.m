function [final_mesh_x, final_mesh_y, omega_x, omega_y] = dynamic_mesh_refinement2D(nx0, x0, x1,ny0, y0, y1, epsilon, function_string)

    count_iter = 0;

    % n0 = initial value for the number of mesh nodes
    computational_domain_x = linspace(x0, x1, nx0);
    computational_domain_x = computational_domain_x';
    computational_domain_y = linspace(y0, y1, ny0);
    computational_domain_y = computational_domain_y';
   
    physical_domain_x = computational_domain_x;    
    physical_domain_y = computational_domain_y;    
    f = str2func(function_string);
    [A,rhs] = f(physical_domain_x, physical_domain_y);
    sol = gmres(A,rhs,50,1e-8,200, diag(diag(A))) ;
    
    discrepancy_x = 1;
    discrepancy_y = 1;
    
    while( discrepancy_x > epsilon || discrepancy_y > epsilon )
        
        %Extraxt strides along x
        omega_x_strides = [];
        
        for stride = 2:ny0-1
            %extract strides of solution along x
            sol_strides_x = sol((stride-1)*nx0+1:stride*nx0);
            gradient_solution = Centered_gradient(physical_domain_x, sol_strides_x);
            omega = monitor_function1D(physical_domain_x, gradient_solution); 
            omega_x_strides = [omega_x_strides; omega'];
        end
        
        temp_physical_domain_x = zeros(size(computational_domain_x,1),1);
        temp_physical_domain_x(1) = x0;
        
        omega_x = [];
        
        for stride = 1:size(omega_x_strides,1)  
            omega_x = [omega_x; max(omega_x_strides(stride,:))]; 
        end
        
        for i = 1:size(omega_x,1)          
            temp_physical_domain_x(i+1) = temp_physical_domain_x(i) + 1/omega_x(i) ;     
        end
        
        temp_physical_domain_x(end) = temp_physical_domain_x(end-1) + 1/omega_x(end-1);
        
        % renormalize 
        temp_physical_domain_x = x0 + temp_physical_domain_x / (temp_physical_domain_x(end) - temp_physical_domain_x(1)) * (x1-x0);
        discrepancy_x = max( physical_domain_x - temp_physical_domain_x ); 
        physical_domain_x = temp_physical_domain_x; 
        
        
        %Extraxt strides along y
        omega_y_strides = [];
        
        for stride = 2:nx0-1
            %extract strides of solution along x
            sol_strides_y = zeros(ny0,1); 
            for index = 1:ny0
                sol_strides_y(index) = sol(stride+(index-1)*nx0);
            end
            gradient_solution = Centered_gradient(physical_domain_y, sol_strides_y);
            omega = monitor_function1D(physical_domain_y, gradient_solution); 
            omega_y_strides = [omega_y_strides; omega'];
        end
        
        temp_physical_domain_y = zeros(size(computational_domain_y,1),1);
        temp_physical_domain_y(1) = y0;
        
        omega_y = [];
        
        for stride = 1:size(omega_y_strides,1)  
            omega_y = [omega_y; max(omega_y_strides(stride,:))]; 
        end
        
        for i = 1:size(omega_y,1)          
            temp_physical_domain_y(i+1) = temp_physical_domain_y(i) + 1/omega_y(i) ;     
        end
        
        temp_physical_domain_y(end) = temp_physical_domain_y(end-1) + 1/omega_y(end-1);
        
        % renormalize 
        temp_physical_domain_y = y0 + temp_physical_domain_y / (temp_physical_domain_y(end) - temp_physical_domain_y(1)) * (y1-y0);
        discrepancy_y = max( physical_domain_y - temp_physical_domain_y ); 
        physical_domain_y = temp_physical_domain_y; 
        
        
        [A,rhs] = f(physical_domain_x, physical_domain_y);
        sol = gmres(A,rhs,50,1e-8,200, diag(diag(A))) ;
        
        display(strcat('Iteration counter: ', num2str(count_iter), ' - discrepancy_x: ', num2str(discrepancy_x), ' - discrepancy_y: ', num2str(discrepancy_y)));
         
    end
    
    final_mesh_x = physical_domain_x; 
    final_mesh_y = physical_domain_y; 

return 
end