function [Per]=redblack(G)

    %reordering of the grid points
    nx=size(G,2);
    ny=size(G,1);
    R=zeros(size(G));

    for j=2:ny-1
        for i=2:nx-1
            if mod(j-1,2)==1 && mod(i-1,2)==1
               R(i,j)=1;
            end
            if mod(j-1,2)==0 && mod(i-1,2)==0
                R(i,j)=1;
            end
        end
    end

    %creation of the permutation matrix
    count=1;
    Per=zeros((n-2)^2,(n-2)^2);
    for j=2:ny-1
        for i=2:nx-1
            if R(i,j)==1
                Per(G(i,j),count)=1;
                count=count+1;
            end
        end
    end

    for j=2:ny-1
        for i=2:nx-1
            if R(i,j)==0
                Per(G(i,j),count)=1;
                count=count+1;
            end
        end
    end

end

