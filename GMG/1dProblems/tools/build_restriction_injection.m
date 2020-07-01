function R = build_restriction_injection( n0, n1, x1, x2 )

    h = (x2 - x1)/(n0 - 1);   
    H = (x2 - x1)/(n1 - 1);
    
    R = sparse(n1,n0);
    
    ratio = H/h;
    
    count = 0;
    
    for i = 1:n1
        R(i, count*ratio + 1) = 1;
        count = count + 1;
    end


return
end