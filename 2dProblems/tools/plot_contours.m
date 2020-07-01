function plot_contours( sol_vector, nx, ny )

    Z = zeros(nx,ny);
    
    for i = 1:nx 
        for j = 1:ny
            Z(i,j) = sol_vector(i + (j-1)*nx);
        end
    end

    [M,c] = contourf(Z, 30);
    c.LineWidth = 1;
    colorbar()
    colormap('jet')

end