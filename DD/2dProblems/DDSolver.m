close all
clear
clc

addpath('./tools')
addpath('../../2dProblems')
addpath('../../2dProblems/tools')

Nx = 101;
Ny = 101;
x1 = 0.0;
x2 = 1.0;
y1 = 0.0;
y2 = 1.0;


ndomains_x = 5;
ndomains_y = 5;

n_domains = ndomains_x * ndomains_y;

overlapping = 0.4;

problem = @(x1,x2,y1,y2,Nx,Ny) fd2d_heat_steady_test02 ( x1,x2,y1,y2,Nx,Ny );

[DD_data, count_contributions] = PrepareDD_Data( Nx, x1, x2, ndomains_x, Ny, y1, y2, ndomains_y, overlapping, problem );    

%Computation of exact solution
[A,rhs] = problem ( x1,x2,y1,y2,Nx,Ny );
exact = A\rhs;

u = zeros(Nx*Ny,1);

error = norm(u-exact)/norm(exact);

count = 0;

while( error > 10^(-6) )
    
    v = zeros(Nx*Ny,1);
    
    for dom = 1 : n_domains
        R = DD_data{dom}.R;
        v( DD_data{dom}.global_indices ) = v( DD_data{dom}.global_indices ) +  DD_data{dom}.A \ (R * ( rhs - A * u) );
    end
    
    u = u + count_contributions.\v;
    
    error = norm(u-exact)/norm(exact);

    display( error )

    count = count + 1;
    
end
