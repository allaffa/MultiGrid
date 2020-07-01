function  [W, DD_data] = colouring(DD_data, global_size, ndx)

    % In 2D the total number of labels is 2 

    ndomains = size(DD_data,1);
    
    W = zeros(global_size,9);
    
    for i = 1:ndomains
        
        grid_row = ceil(i/ndx);
        ypos = mod(grid_row,3);
        
        if( mod(i,ndx)~=0 )
            xpos = mod(i - floor(i/ndx) * ndx,3);
        else
            xpos = xpos = mod(i - ((ypos - 1)*ndx)), 3);
        end
        
        if (xpos > 0) && (ypos > 0)        
            label = (ypos - 1) * 3 + xpos;             
        elseif (xpos > 0) && (ypos == 0)   
            label = 6 + xpos; 
        elseif (xpos == 0) && (ypos > 0)  
            label = ypos * 3;
        else 
            label = 9;
        end
        
        W(DD_data{i}.global_indices,label) = ones(size(DD_data{i}.global_indices,1),1);
        DD_data{i}.colouring = label;
        
    end

end
