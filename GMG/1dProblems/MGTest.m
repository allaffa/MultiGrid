close all
clear
clc

addpath('./tools')
addpath('../../1dProblems')
addpath('../../1dProblems/tools')

N = 101;
x1 = 0.0;
x2 = 1.0;

Hratio = 10;

relaxing = 'Jacobi';
omega = 1;  

problem = @(N,x1,x2 ) fd1d_bvp_test01 ( N,x1,x2 );

MGData = PrepareMGData( N, x1,x2, Hratio, relaxing, problem );   

%Computation of exact solution
[A,rhs] = problem ( N,x1,x2 );
exact = A\rhs;

num_cycles = 70;
error = zeros(num_cycles,1);
trend = zeros(num_cycles,1);

x0 = zeros(size(A,1),1);

plot_theoretical_curve = 0;

if strcmp(relaxing, 'Jacobi')
    
    [omega, ratio] = compute_convergence_ratio( Hratio, N );
    plot_theoretical_curve = 1;
    
elseif strcmp(relaxing, 'GaussSeidel')
    
    ratio = 1/sqrt(5);
    plot_theoretical_curve = 1;
    
end

for cycle = 1:num_cycles
    
    [u] = V_cycle( 1, omega, relaxing, MGData, x0 );
    error(cycle) = norm( exact - u )/norm( exact );
    
    if plot_theoretical_curve == 1;
        trend(cycle) = (ratio) ^ (cycle-1);
    end
        
    x0 = u;
    
end

figure()
semilogy( error, 'r', 'linewidth',3 )
hold on

if plot_theoretical_curve == 1;
    semilogy( trend, '--', 'linewidth',3 )
end
    
set(gca, 'fontsize', 18)
legend('MG error', 'Theoretical error descent')

