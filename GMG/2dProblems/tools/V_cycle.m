function [u] = V_cycle( nu1, omega, relaxing, MGData, x0 )
    
    MGData{1}.u = x0;
    
    if(strcmp(relaxing, 'Jacobi'))
    
        % Relaxation at the finest level
        for relax = 1:nu1
            z = MGData{1}.A * MGData{1}.u;
            MGData{1}.u = (1-omega) * MGData{1}.u  +  ...
                omega * ( MGData{1}.u - MGData{1}.Relax \ z + MGData{1}.Relax \ MGData{1}.b );

        end   
    
    else
        
        % Relaxation at the finest level
        for relax = 1:nu1
            z = MGData{1}.A * MGData{1}.u;
            MGData{1}.u = MGData{1}.u - MGData{1}.Relax \ z + MGData{1}.Relax \ MGData{1}.b ;

        end 
        
    end

    %Computation of the residual
    res = MGData{1}.b - MGData{1}.A * MGData{1}.u;
    
    %Restriciton of the residual to the coarse level
    MGData{2}.b = MGData{1}.R * res;

    % Direct solving the residual equation at the coarse level
    MGData{2}.u = MGData{2}.A\MGData{2}.b;

    %Prolungation of the error to the fine level
    update = MGData{2}.P * MGData{2}.u;
    
    %Solution updating
    u = MGData{1}.u + update;
    MGData{1}.u = u;
             
return 
end