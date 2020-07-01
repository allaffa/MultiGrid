close all
clear
clc

addpath('./tools')
addpath('../../1dProblems')
addpath('../../1dProblems/tools')

N = 101;
x1 = 0.0;
x2 = 1.0;

n_domains = 2;

% overlapping is expressed in percentage with respect to the size of the
% single subdomain. Therefore it must be a number between 0 and 1

overlapping = 0.2;

problem = @(N,x1,x2 ) fd1d_bvp_test01 ( N,x1,x2 );

[DD_data, count_contributions] = PrepareDD_Data( N, x1,x2, n_domains, overlapping, problem );   

%Computation of exact solution
[A,rhs] = problem ( N,x1,x2 );
exact = A\rhs;

u = zeros(N,1);

error = norm(u-exact)/norm(exact);

count = 0;

while( error > 10^(-6) )
    
    v = zeros(N,1);
    
    for dom = 1 : n_domains
        R = DD_data{dom}.R;
        v( DD_data{dom}.global_indices ) = v( DD_data{dom}.global_indices ) +  DD_data{dom}.A \ (R * ( rhs - A * u) );
    end
    
    u = u + count_contributions.\v;
    
    error = norm(u-exact)/norm(exact);

    display( error )

    count = count + 1;
    
end
