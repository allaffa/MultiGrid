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

%% Creation of covariate matrix 

NDOMAINS = [ 2 3 4 5 6 7 8 9 10 11 ]';

ratiox = 50;
ratioy = 50;

Z = zeros(length(NDOMAINS),4);

TIME_dd = zeros(size(NDOMAINS));
TIME_p = zeros(size(NDOMAINS));
TIME_ac = zeros(size(NDOMAINS));

TIME_dd_modified2 = zeros(size(NDOMAINS));
TIME_p_modified2 = zeros(size(NDOMAINS));
TIME_pbar_modified2 = zeros(size(NDOMAINS));
TIME_ac_modified2 = zeros(size(NDOMAINS));
TIME_ibar_modified2 = zeros(size(NDOMAINS));

for i = 1:length(NDOMAINS)
    
    Z(i,1) = 1;
    Z(i,2) = Z(i,1) * (NDOMAINS(i) * ratiox * NDOMAINS(i) * ratioy);
    Z(i,3) = Z(i,2) * (NDOMAINS(i) * ratiox * NDOMAINS(i) * ratioy);
    Z(i,4) = Z(i,3) * (NDOMAINS(i) * ratiox * NDOMAINS(i) * ratioy);
  
end


%% Definition of global variables to measure performance

global DD_time Ac_time P_time ...
    DD_time_modified2 Ac_time_modified2 P_time_modified2 Pbar_time_modified2 Ibar_time_modified2

for k = 1:length(NDOMAINS)

    %% Domain setting

    Nx = NDOMAINS(k) * ratiox + 1;
    Ny = NDOMAINS(k) * ratioy + 1;
    x1 = 0.0;
    x2 = 5.0;
    y1 = 0.0;
    y2 = 5.0;

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

    %% Krylov method setting

    tol = 10^(-8);
    maxit = 10^(3);
    params(1) = tol; params(2) = maxit;  % tolerance and max iterations

    %% Number of runs to averagre over
    
    number_runs = 100;

%% Sequel of subdomain sizes to use

    display( strcat('Number of subdomains: ',  num2str(NDOMAINS(k)*NDOMAINS(k)) ) ) ;
        
    DD_time = 0.0;
    Ac_time = 0.0;
    P_time = 0.0;
    
    DD_time_modified2 = 0.0;
    Ac_time_modified2 = 0.0;
    P_time_modified2 = 0.0;
    Pbar_time_modified2 = 0.0;
    Ibar_time_modified2 = 0.0;
    
    MGData = PrepareMGData( Nx, x1,x2, ratiox, Ny, y1, y2, ratioy, relaxing, problem );   

    ndomains_x = NDOMAINS(k);
    ndomains_y = NDOMAINS(k);
    
    n_domains = ndomains_x * ndomains_y;

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
    
    TIME_dd(k) = (DD_time/n_domains)/number_runs;
    TIME_p(k) = (P_time/n_domains)/number_runs;
    TIME_ac(k) = (Ac_time)/number_runs;
         
    TIME_dd_modified2(k) = (DD_time_modified2/n_domains)/number_runs;
    TIME_p_modified2(k) = (P_time_modified2/n_domains)/number_runs;
    TIME_pbar_modified2(k) = (Pbar_time_modified2/n_domains)/number_runs;
    TIME_ac_modified2(k) = (Ac_time_modified2)/number_runs;
    TIME_ibar_modified2(k) = (Ibar_time_modified2/n_domains)/number_runs;
    


end

coeffs_dd = (Z' * Z)\(Z' * TIME_dd);
coeffs_p = (Z' * Z)\(Z' * TIME_p);
coeffs_ac = (Z' * Z)\(Z' * TIME_ac);

coeffs_dd_modified2 = (Z' * Z)\(Z' * TIME_dd_modified2);
coeffs_p_modified2 = (Z' * Z)\(Z' * TIME_p_modified2);
coeffs_pbar_modified2 = (Z' * Z)\(Z' * TIME_pbar_modified2);
coeffs_ac_modified2 = (Z' * Z)\(Z' * TIME_ac_modified2);
coeffs_ibar_modified2 = (Z' * Z)\(Z' * TIME_ibar_modified2);


reg_coeffs_dd = @(x) coeffs_dd(1) + coeffs_dd(2) * x + coeffs_dd(3) * x .* x + coeffs_dd(4) * x .* x .* x;
reg_coeffs_p = @(x) coeffs_p(1) + coeffs_p(2) * x + coeffs_p(3) * x .* x + coeffs_p(4) * x .* x .* x;
reg_coeffs_ac = @(x) coeffs_ac(1) + coeffs_ac(2) * x + coeffs_ac(3) * x .* x + coeffs_ac(4) * x .* x .* x;

reg_coeffs_dd_modified2 = @(x) coeffs_dd_modified2(1) + coeffs_dd_modified2(2) * x + coeffs_dd_modified2(3) * x .* x + coeffs_dd_modified2(4) * x .* x .* x;
reg_coeffs_p_modified2 = @(x) coeffs_p_modified2(1) + coeffs_p_modified2(2) * x + coeffs_p_modified2(3) * x .* x + coeffs_p_modified2(4) * x .* x .* x;
reg_coeffs_pbar_modified2 = @(x) coeffs_pbar_modified2(1) + coeffs_pbar_modified2(2) * x + coeffs_pbar_modified2(3) * x .* x + coeffs_pbar_modified2(4) * x .* x .* x;
reg_coeffs_ac_modified2 = @(x) coeffs_ac_modified2(1) + coeffs_ac_modified2(2) * x + coeffs_ac_modified2(3) * x .* x + coeffs_ac_modified2(4) * x .* x .* x;
reg_coeffs_ibar_modified2 = @(x) coeffs_ibar_modified2(1) + coeffs_ibar_modified2(2) * x + coeffs_ibar_modified2(3) * x .* x + coeffs_ibar_modified2(4) * x .* x .* x;

figure()
semilogy((NDOMAINS * ratiox .* NDOMAINS * ratioy), reg_coeffs_dd((NDOMAINS * ratiox .* NDOMAINS * ratioy)), 'linewidth',5);
hold on
semilogy((NDOMAINS * ratiox .* NDOMAINS * ratioy), TIME_dd, 'o', 'linewidth', 10);
semilogy((NDOMAINS* ratiox .* NDOMAINS * ratioy), reg_coeffs_p((NDOMAINS* ratiox .* NDOMAINS * ratioy)), 'linewidth',5);
semilogy((NDOMAINS* ratiox .* NDOMAINS * ratioy), TIME_p,'o', 'linewidth', 10);
semilogy((NDOMAINS* ratiox .* NDOMAINS * ratioy), reg_coeffs_ac((NDOMAINS* ratiox .* NDOMAINS * ratioy)), 'linewidth',5);
semilogy((NDOMAINS* ratiox .* NDOMAINS * ratioy), TIME_ac,'o', 'linewidth', 10);
set(gca, 'fontsize', 25)
title('Timing for Standard approach')
xlabel('Problem size')
ylabel('CPU Time (s)')
legend('DD - regression', 'DD - sample points', 'Prolongation and restriciton - sample points', 'Prolongation and restriction - regression', 'Coarse grid solver - regression', ...
    'Coarse grid solver - sample points')

figure()
semilogy((NDOMAINS * ratiox .* NDOMAINS * ratioy), reg_coeffs_dd_modified2((NDOMAINS * ratiox .* NDOMAINS * ratioy)), 'linewidth',5);
hold on
semilogy((NDOMAINS * ratiox .* NDOMAINS * ratioy), TIME_dd_modified2, 'o', 'linewidth', 10);
semilogy((NDOMAINS* ratiox .* NDOMAINS * ratioy), reg_coeffs_p_modified2((NDOMAINS* ratiox .* NDOMAINS * ratioy)), 'linewidth',5);
semilogy((NDOMAINS* ratiox .* NDOMAINS * ratioy), TIME_p_modified2,'o', 'linewidth', 10);
semilogy((NDOMAINS* ratiox .* NDOMAINS * ratioy), reg_coeffs_pbar_modified2((NDOMAINS* ratiox .* NDOMAINS * ratioy)), 'linewidth',5);
semilogy((NDOMAINS* ratiox .* NDOMAINS * ratioy), TIME_pbar_modified2,'o', 'linewidth', 10);
semilogy((NDOMAINS* ratiox .* NDOMAINS * ratioy), reg_coeffs_ac_modified2((NDOMAINS* ratiox .* NDOMAINS * ratioy)), 'linewidth',5);
semilogy((NDOMAINS* ratiox .* NDOMAINS * ratioy), TIME_ac_modified2,'o', 'linewidth', 10);
set(gca, 'fontsize', 25)
title('Timing for Modified approach')
xlabel('Problem size')
ylabel('CPU Time (s)')
legend('DD - regression', 'DD - sample points', 'Restriciton - sample points', 'Restriction - regression', 'Smoother prolongation - regression',...
    'Smoother prolongation - sample points', 'Coarse grid solver - regression', 'Coarse grid solver - sample points')



save('Regression2.mat');
