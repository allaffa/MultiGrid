close all
clear
%clc

%% Inclusion of the paths to use the tools from GMG and DD

addpath('../../StartUpSetting')
addpath('../../2dProblems')
addpath('../../2dProblems/tools')
addpath('../../GMG/2dProblems/tools')
addpath('../../DD/2dProblems/tools')
addpath('../../LinearSolvers')
addpath('../Preconditioners')

%% Set up multicore policy

startup


%% Domain setting

Nx = 301;
Ny = 301;
x1 = 0.0;
x2 = 3.0;
y1 = 0.0;
y2 = 3.0;

%% Problem definition

display('Construction of the problem on the fine mesh')

problem = @(x1,x2,y1,y2,Nx,Ny) fd2d_heat_steady_test02 ( x1,x2,y1,y2,Nx,Ny );

%Computation of exact solution
[A,rhs] = problem ( x1,x2,y1,y2,Nx,Ny );
exact = A\rhs;

%% Mutligrid relaxation
relaxing = 'Jacobi';

%% Domain decomposition overlapping 
overlapping = 0.0;

%% Definition of global variables to measure performance

global DD_time Ac_time P_time ...
    DD_time_modified2 Ac_time_modified2 P_time_modified2 Pbar_time_modified2 Ibar_time_modified2

%% Krylov method setting

tol = 10^(-8);
maxit = 10^(3);
params(1) = tol; params(2) = maxit;  % tolerance and max iterations

%% Number of runs to averagre over
number_runs = 100;

%% Sequel of subdomain sizes to use

%SIZE = [ 2 4 8 16 20 25 40 50 80 100 200 ]';
SIZE = [ 3 6 10 20 30 50 60 100 150 ]';
NDOMAINS = zeros( size(SIZE) );
Z = zeros(length(SIZE),4);

TIME_dd = zeros(size(SIZE));
TIME_p = zeros(size(SIZE));
TIME_ac = zeros(size(SIZE));

TIME_dd_modified2 = zeros(size(SIZE));
TIME_p_modified2 = zeros(size(SIZE));
TIME_pbar_modified2 = zeros(size(SIZE));
TIME_ac_modified2 = zeros(size(SIZE));
TIME_ibar_modified2 = zeros(size(SIZE));

for i = 1:length(SIZE)
    
    Z(i,1) = 1;
    Z(i,2) = Z(i,1) * SIZE(i);
    Z(i,3) = Z(i,2) * SIZE(i);
    Z(i,4) = Z(i,3) * SIZE(i);
   
end


for i = 1:length(SIZE)
    
    display( strcat('Size of the domains: ', num2str(SIZE(i) )) )
        
    DD_time = 0.0;
    Ac_time = 0.0;
    P_time = 0.0;
    
    DD_time_modified2 = 0.0;
    Ac_time_modified2 = 0.0;
    P_time_modified2 = 0.0;
    Pbar_time_modified2 = 0.0;
    Ibar_time_modified2 = 0.0;
    
    ratiox = SIZE(i);
    ratioy = SIZE(i);

    MGData = PrepareMGData( Nx, x1,x2, ratiox, Ny, y1, y2, ratioy, relaxing, problem );   

    ndomains_x = (Nx-1)/ratiox;
    ndomains_y = (Ny-1)/ratioy;
    
    n_domains = ndomains_x * ndomains_y;
    NDOMAINS(i) = n_domains;

    [DD_data, count_contributions] = PrepareDD_Data( Nx, x1, x2, ndomains_x, Ny, y1, y2, ndomains_y, overlapping, problem );  

    colour = 0;
    start = cputime;
    [Smoothed_Prolungation] = Hybrid2Schwarz_Smoothed_interpolation(A, MGData, DD_data, count_contributions, colour);
    finish = cputime;
    Ibar_time_modified2 = Ibar_time_modified2 + (finish - start);

    P = @(x) Hybrid2Schwarz_Preconditioner(A, MGData, DD_data, count_contributions, x);
    Pmodified2 = @(x) Hybrid2Schwarz_Preconditioner_modified2(A, Smoothed_Prolungation, MGData, DD_data, count_contributions, x, colour);
    
    for runs = 1:number_runs
    
        initialguess = rand(size(A,1),1);

        [u_prec,flag_prec,iter_prec,absres_prec] = mygmres(A,P,rhs,params,initialguess);
        [u_prec_modified2,flag_prec_modified2,iter_prec_modified2,absres_prec_modified2] = mygmres(A,Pmodified2,rhs,params,initialguess);

    end
    
    TIME_dd(i) = (DD_time/n_domains)/number_runs;
    TIME_p(i) = (P_time/n_domains)/number_runs;
    TIME_ac(i) = (Ac_time)/number_runs;
         
    TIME_dd_modified2(i) = (DD_time_modified2/n_domains)/number_runs;
    TIME_p_modified2(i) = (P_time_modified2/n_domains)/number_runs;
    TIME_pbar_modified2(i) = (Pbar_time_modified2/n_domains)/number_runs;
    TIME_ac_modified2(i) = (Ac_time_modified2)/number_runs;
    TIME_ibar_modified2(i) = (Ibar_time_modified2/n_domains)/number_runs;
    
end

coeffs_dd = (Z' * Z)\(Z' * TIME_dd);
coeffs_p = (Z' * Z)\(Z' * TIME_p);
coeffs_ac = (Z' * Z)\(Z' * TIME_ac);

coeffs_dd_modified2 = (Z' * Z)\(Z' * TIME_dd_modified2);
coeffs_p_modified2 = (Z' * Z)\(Z' * TIME_p_modified2);
coeffs_pbar_modified2 = (Z' * Z)\(Z' * TIME_pbar_modified2);
coeffs_ac_modified2 = (Z' * Z)\(Z' * TIME_ac_modified2);
coeffs_ibar_modified2 = (Z' * Z)\(Z' * TIME_ibar_modified2);


reg_coeffs_dd = @(x) coeffs_dd(1) + coeffs_dd(2) * x + coeffs_dd(3) * x.* x + coeffs_dd(4) * x.* x .* x;
reg_coeffs_p = @(x) coeffs_p(1) + coeffs_p(2) * x + coeffs_p(3) * x.* x + coeffs_p(4) * x .* x.* x;
reg_coeffs_ac = @(x) coeffs_ac(1) + coeffs_ac(2) * x + coeffs_ac(3) * x .* x + coeffs_ac(4) * x.* x.* x;

reg_coeffs_dd_modified2 = @(x) coeffs_dd_modified2(1) + coeffs_dd_modified2(2) * x + coeffs_dd_modified2(3) * x .* x + coeffs_dd_modified2(4) * x .* x .* x;
reg_coeffs_p_modified2 = @(x) coeffs_p_modified2(1) + coeffs_p_modified2(2) * x + coeffs_p_modified2(3) * x .* x + coeffs_p_modified2(4) * x .* x .* x;
reg_coeffs_pbar_modified2 = @(x) coeffs_pbar_modified2(1) + coeffs_pbar_modified2(2) * x + coeffs_pbar_modified2(3) * x .* x + coeffs_pbar_modified2(4) * x .* x .* x;
reg_coeffs_ac_modified2 = @(x) coeffs_ac_modified2(1) + coeffs_ac_modified2(2) * x + coeffs_ac_modified2(3) * x .* x + coeffs_ac_modified2(4) * x .* x .* x;
reg_coeffs_ibar_modified2 = @(x) coeffs_ibar_modified2(1) + coeffs_ibar_modified2(2) * x + coeffs_ibar_modified2(3) * x .* x + coeffs_ibar_modified2(4) * x .* x .* x;

figure()
plot(SIZE, reg_coeffs_dd(SIZE), 'linewidth',5);
hold on
plot(SIZE, TIME_dd, 'o', 'linewidth', 10);
plot(SIZE, reg_coeffs_dd_modified2(SIZE), '--', 'linewidth', 5);
plot(SIZE, TIME_dd_modified2, 'o', 'linewidth', 10);
set(gca, 'fontsize', 25)
title('Subdomain DD solver')
xlabel('Subdomain size')
ylabel('CPU Time (s)')
legend('Standard approach - regression', 'Standard approach - sample points', 'Modified approach 2 - regression', 'Modified approach 2 - sample points')

figure()
plot(SIZE, reg_coeffs_p(SIZE), 'linewidth',5);
hold on
plot(SIZE, TIME_p,'o', 'linewidth', 10);
plot(SIZE, reg_coeffs_p_modified2(SIZE), '--', 'linewidth', 5);
plot(SIZE, TIME_p_modified2, 'o', 'linewidth', 10);
set(gca, 'fontsize', 25)
title('Interpolation and restriction')
xlabel('Subdomain size')
ylabel('CPU Time (s)')
legend('Standard approach - regression', 'Standard approach - sample points', 'Modified approach 2 - regression', 'Modified approach 2 - sample points')


figure()
plot(SIZE, reg_coeffs_ac(SIZE), 'linewidth',5);
hold on
plot(SIZE, TIME_ac, 'o', 'linewidth', 10);
plot(SIZE, reg_coeffs_ac_modified2(SIZE), '--', 'linewidth', 5);
plot(SIZE, TIME_ac_modified2, 'o', 'linewidth', 10);
set(gca, 'fontsize', 25)
title('Multigrid solver')
xlabel('Subdomain size')
ylabel('CPU Time (s)')
legend('Standard approach - regression', 'Standard approach - sample points', 'Modified approach 2 - regression', 'Modified approach 2 - sample points')

figure()
plot(SIZE, reg_coeffs_pbar_modified2(SIZE), 'linewidth',5);
hold on
plot(SIZE, TIME_pbar_modified2, 'o', 'linewidth', 10);
set(gca, 'fontsize', 25)
title('Smoothed interpolation operator (mat-vec)')
xlabel('Subdomain size')
ylabel('CPU Time (s)')
legend('Modified approach 2 - regression', 'Modified approach 2 - sample points')

figure()
plot(SIZE, reg_coeffs_ibar_modified2(SIZE), 'linewidth',5);
hold on
plot(SIZE, TIME_ibar_modified2, 'o', 'linewidth', 10);
set(gca, 'fontsize', 25)
title('Smoothed interpolation operator (construction)')
xlabel('Subdomain size')
ylabel('CPU Time (s)')
legend('Modified approach 2 - regression', 'Modified approach 2 - sample points')


save('Regression.mat');
