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

%% For the construction of the MG and DD data structures we do't need to measure the performance, so I use all the cores possible for this

maxNumCompThreads('automatic')

%% Domain setting

Nx = 301;
Ny = 301;
x1 = 0.0;
x2 = 3.0;
y1 = 0.0;
y2 = 3.0;

%% Problem definition

display('Construction of the problem on the fine mesh')

problem = @(x1,x2,y1,y2,Nx,Ny) fd2d_heat_steady_test04 ( x1,x2,y1,y2,Nx,Ny );

%Computation of exact solution
[A,rhs] = problem ( x1,x2,y1,y2,Nx,Ny );
exact = A\rhs;

RATIO = [2 4 5 6 10 15 20 30 50 60 100]';

nnz_A    = zeros(size(RATIO));
nnz_L    = zeros(size(RATIO));
nnz_U    = zeros(size(RATIO));
nnz_Ibar = zeros(size(RATIO));
nnz_I = zeros(size(RATIO));

for settings = 1:length(RATIO)

    % level of coarsening for the coarse mesh
    ratiox = RATIO(settings);
    ratioy = RATIO(settings);

    % Relaxation technique used
    relaxing = 'Jacobi';

    display('Construction of the MG data structure')

    % In the data structure, 1 is the index for the fine level, 2 for the
    % coarse level
    MGData = PrepareMGData( Nx, x1,x2, ratiox, Ny, y1, y2, ratioy, relaxing, problem );   


    % Setting the number of domains
    ndomains_x = (Nx - 1) / ratiox;
    ndomains_y = (Ny - 1) / ratioy;

    n_domains = ndomains_x * ndomains_y;

    % overlapping is expressed in percentage with respect to the size of the
    % single subdomain. Therefore it must be a number between 0 and 1
    overlapping = 0.0;

    display('Construction of the DD data structure')

    [DD_data, count_contributions] = PrepareDD_Data( Nx, x1, x2, ndomains_x, Ny, y1, y2, ndomains_y, overlapping, problem );  

    display( 'Construction of the smoothed preconditioner' )

    start = cputime;
    [Smoothed_Prolungation] = Hybrid2Schwarz_Smoothed_interpolation(A, MGData, DD_data, count_contributions);
    finish = cputime;
    
    nnz_A(settings)    = nnz( MGData{2}.A );
    nnz_L(settings)    = nnz( MGData{2}.L );
    nnz_U(settings)    = nnz( MGData{2}.U );
    nnz_I(settings)    = nnz( MGData{2}.P );
    nnz_Ibar(settings) = nnz( Smoothed_Prolungation );
    
end


figure()
semilogy(RATIO, nnz_L/nnz_A, 'b-o', 'linewidth', 3);
hold on
semilogy(RATIO, nnz_U/nnz_A, 'g-o', 'linewidth', 3);
semilogy(RATIO, nnz_Ibar/nnz_I, 'r-o', 'linewidth', 3);
set(gca, 'fontsize', 25);
title('Relative Fill-in')
xlabel('Coarsening ratio H/h')
ylabel('Fill-in')
%legend('nnz(L_H)/nnz(A_H)', 'nnz(U_H)/nnz(A_H)', 'nnz(Ibar)/nnz(I)')


save('sparsity_check.mat')