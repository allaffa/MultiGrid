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

ratiox = 2;
ratioy = 10;

relaxing = 'Jacobi';
omega = 1;

problem = @(x1,x2,y1,y2,Nx,Ny) fd2d_heat_steady_test01 ( x1,x2,y1,y2,Nx,Ny );

MGData = PrepareMGData( Nx, x1,x2, ratiox, Ny, y1, y2, ratioy, relaxing, problem );   

%Computation of exact solution
[A,rhs] = problem ( x1,x2,y1,y2,Nx,Ny );
exact = A\rhs;

num_cycles = 70;
error = zeros(num_cycles,1);
trend = zeros(num_cycles,1);

x0 = zeros(size(A,1),1);

plot_theoretical_curve = 0;

if strcmp(relaxing, 'Jacobi')
    
    [omega, ratio] = compute_convergence_ratio( ratiox, Nx, ratioy, Ny );
    plot_theoretical_curve = 1;
    
elseif strcmp(relaxing, 'GaussSeidel') && ratiox == 2 && ratioy == 2
    
    ratio = 1/2;
    plot_theoretical_curve = 1;
    
end

for cycle = 1:num_cycles
    
    [u] = V_cycle( 1, omega, relaxing, MGData, x0 );
    error(cycle) = norm( exact - u )/norm( exact );
    
    if plot_theoretical_curve == 1
        trend(cycle) = (ratio) ^ (cycle-1);
    end
    
    x0 = u;
    
end

figure()
semilogy( error, 'r', 'linewidth',3 )
hold on

if plot_theoretical_curve == 1
    semilogy( trend, '--', 'linewidth',3 )
end

set(gca, 'fontsize', 18)
legend('MG error', 'Theoretical error descent')

