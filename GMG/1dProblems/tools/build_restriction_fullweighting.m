function R = build_restriction_fullweighting( n0, n1, x1, x2 )

    h = (x2 - x1)/(n0 - 1);   
    H = (x2 - x1)/(n1 - 1);
    
    R = sparse(n1,n0);
    
    ratio = H/h;
    
    denominator = ratio * ratio ;   
    
    
    % The extension of a stencil for the full weighting reduction is 
    %
    %   ratio + 2 * (ratio-1)*ratio/2  (Summation Gauss formula)
    %
    
    
    for i = 2:n1-1
        
        center = (i-1)*ratio + 1;
        
        R( i, center ) = ratio;
        
        for j = 1:(ratio-1)
            R( i, center-j ) = ratio-j;
        end
        
        for j = 1:(ratio-1)
            R( i, center+j ) = ratio-j;
        end
        
    end

    R(1,1) = denominator;
    R(n1,n0) = denominator;
    
    R = 1/denominator * R;

return
end