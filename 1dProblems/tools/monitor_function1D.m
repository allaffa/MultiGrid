function monitor_values = monitor_function1D(xh, gradient)

    assert(size(xh,1)==(size(gradient,1))+2);
    monitor_values = zeros(size(gradient,1),1);
    
    % We only compute the gradient for internal nodes. 
    % In case we have to refine, we will refine also on the adjacent area 
    % on the boundary

    for i = 1:size(gradient,1)
        
        monitor_values(i) = sqrt( 1+gradient(i)^2 );
        
    end

    return;
    
end

