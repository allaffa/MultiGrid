close all
clear
clc

%% Inclusion of the paths to use the tools from GMG and DD

addpath('../../2dProblems')
addpath('../../GMG/2dProblems/tools')
addpath('../../DD/2dProblems/tools')

%% Domain setting

Nx = 101;
Ny = 101;
x1 = 0.0;
x2 = 1.0;
y1 = 0.0;
y2 = 1.0;

%% Problem definition

problem = @(x1,x2,y1,y2,Nx,Ny) fd2d_heat_steady_test01 ( x1,x2,y1,y2,Nx,Ny );

%Computation of exact solution
[A,rhs] = problem ( x1,x2,y1,y2,Nx,Ny );
exact = A\rhs;

%% Multigrid setting

% level of coarsening for the coarse mesh
ratiox = 10;
ratioy = 10;

% Relaxation technique used
relaxing = 'Jacobi';
omega = 1;  

% In the data structure, 1 is the index for the fine level, 2 for the
% coarse level
MGData = PrepareMGData( Nx, x1,x2, ratiox, Ny, y1, y2, ratioy, relaxing, problem );   

%% Domain Decompodition setting

% Setting the number of domains
ndomains_x = 20;
ndomains_y = 20;

n_domains = ndomains_x * ndomains_y;

% overlapping is expressed in percentage with respect to the size of the
% single subdomain. Therefore it must be a number between 0 and 1
overlapping = 0.0;

[DD_data, count_contributions] = PrepareDD_Data( Nx, x1, x2, ndomains_x, Ny, y1, y2, ndomains_y, overlapping, problem );  

%% Two Level Schwarz solver

% Initial guess
u = zeros(size(A,1),1);

% tolerance for the error check
tol = 10^(-6);

error = norm(u-exact)/norm(exact);

count = 0;

while( error > 10^(-6) )
    
    %computation of the residual
    resF = rhs - A * u;
    
    %Restriciton of the residual to the coarse level
    resC = MGData{1}.R * resF;

    % Direct solving the residual equation at the coarse level
    updateC = MGData{2}.A \ resC;

    %Prolungation of the error to the fine level
    updateF = MGData{2}.P * updateC;
    
    %Solution updating
    u = u + updateF;

    v = zeros(size(A,1),1);
    
    for dom = 1 : n_domains
        R = DD_data{dom}.R;
        v( DD_data{dom}.global_indices ) = v( DD_data{dom}.global_indices ) +  DD_data{dom}.A \ (R * ( rhs - A * u) );
    end
    
    u = u + count_contributions.\v;
    
    error = norm(u-exact)/norm(exact);

    display( error )

    count = count + 1;
    
end



