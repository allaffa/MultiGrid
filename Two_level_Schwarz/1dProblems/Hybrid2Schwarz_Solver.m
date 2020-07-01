close all
clear
clc

%% Inclusion of the paths to use the tools from GMG and DD

addpath('../../1dProblems')
addpath('../../GMG/1dProblems/tools')
addpath('../../DD/1dProblems/tools')

%% Domain setting

N = 101;
x1 = 0.0;
x2 = 1.0;

%% Problem definition

problem = @(N,x1,x2 ) fd1d_bvp_test01 ( N,x1,x2 );

%Computation of exact solution
[A,rhs] = problem ( N,x1,x2 );
exact = A\rhs;

%% Multigrid setting

% level of coarsening for the coarse mesh
Hratio = 10;

% Relaxation technique used
relaxing = 'Jacobi';
omega = 1;  

% In the data structure, 1 is the index for the fine level, 2 for the
% coarse level
MGData = PrepareMGData( N, x1,x2, Hratio, relaxing, problem );   

%% Domain Decompodition setting

% Setting the number of domains
n_domains = 2;

% overlapping is expressed in percentage with respect to the size of the
% single subdomain. Therefore it must be a number between 0 and 1
overlapping = 0.2;

[DD_data, count_contributions] = PrepareDD_Data( N, x1,x2, n_domains, overlapping, problem );   

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



