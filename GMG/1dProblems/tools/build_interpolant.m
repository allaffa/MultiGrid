function P = build_interpolant( n0, n1, x1, x2 )

    h = (x2 - x1)/(n0 - 1);   
    H = (x2 - x1)/(n1 - 1);
    
    P = sparse(n0,n1);
    
    ratio = H/h;
    
    for j = 1 : n1-1
                
        for m = 0:ratio-1
            
            P( (j - 1) * ratio + m + 1, j )   = ratio - m;
            P( (j - 1) * ratio + m + 1, j+1 ) = m;
            
        end
                
    end
    
    P( n0,n1 ) = ratio;
    
    P = 1/ratio * P;
    
return
end