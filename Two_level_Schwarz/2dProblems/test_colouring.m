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

Nx = 101;
Ny = 101;
x1 = 0.0;
x2 = 4.0;
y1 = 0.0;
y2 = 4.0;

%% Problem definition

display('Construction of the problem on the fine mesh')

problem = @(x1,x2,y1,y2,Nx,Ny) fd2d_heat_steady_test04 ( x1,x2,y1,y2,Nx,Ny );

%Computation of exact solution
[A,rhs] = problem ( x1,x2,y1,y2,Nx,Ny );
exact = A\rhs;

%% Multigrid setting

% level of coarsening for the coarse mesh
ratiox = 10;
ratioy = 10;

% Relaxation technique used
relaxing = 'Jacobi';

display('Construction of the MG data structure')

% In the data structure, 1 is the index for the fine level, 2 for the
% coarse level
MGData = PrepareMGData( Nx, x1,x2, ratiox, Ny, y1, y2, ratioy, relaxing, problem );   

%% Domain Decomposition setting

% Setting the number of domains
ndomains_x = 10;
ndomains_y = 10;

n_domains = ndomains_x * ndomains_y;

% overlapping is expressed in percentage with respect to the size of the
% single subdomain. Therefore it must be a number between 0 and 1
overlapping = 0.0;

display('Construction of the DD data structure')

[DD_data, count_contributions] = PrepareDD_Data( Nx, x1, x2, ndomains_x, Ny, y1, y2, ndomains_y, overlapping, problem );  

%% Set up multicore policy

startup

%% Explicit construction of the smoothed interpolation

display( 'Construction of the smoothed preconditioner' )

colour = 0;

start = cputime;
[Smoothed_Prolungation_nocolor] = Hybrid2Schwarz_Smoothed_interpolation(A, MGData, DD_data, count_contributions, colour);
finish = cputime;

Ibar_time = finish - start;

% Vassilevskj approach with explicitly built smoothed inteprolator
Pmodified2 = @(x) Hybrid2Schwarz_Preconditioner_modified2(A, Smoothed_Prolungation_nocolor, MGData, DD_data, count_contributions, x, colour);

colour = 1;

[W, DD_data] = colouring(DD_data, Nx*Ny, ndomains_x);

[Smoothed_Prolungation_color] = Hybrid2Schwarz_Smoothed_interpolation(A, MGData, DD_data, count_contributions, colour);

% Vassilevskj approach with explicitly built smoothed inteprolator using
% colouring
Pmodified2_colouring = @(x) Hybrid2Schwarz_Preconditioner_modified2(A, Smoothed_Prolungation_color, MGData, DD_data, count_contributions, x, colour);

%% Call to Krylov solver 

display('Solver called');
number_runs = 100;

tot_prec_modified2_nocolouring = 0.0;
tot_prec_modified2_colouring = 0.0;

tol = 10^(-8);
%restart = 10; Since we are using Full GMRES actually we do not need the
%restart parameter 
maxit = 10^(3);

params(1) = tol; params(2) = maxit;  % tolerance and max iterations


for runs = 1:number_runs
    
    initialguess = rand(size(A,1),1);

    start_prec_modified2_nocolouring = cputime;
    [u_prec_modified2_nocolouring,flag_prec_modified2_nocolouring,iter_prec_modified2_nocolouring,absres_prec_modified2_nocolouring] = mygmres(A,Pmodified2,rhs,params,initialguess);
    finish_prec_modified2_nocolouring = cputime;

    tot_prec_modified2_nocolouring = tot_prec_modified2_nocolouring + (finish_prec_modified2_nocolouring - start_prec_modified2_nocolouring);
    
     start_prec_modified2_colouring = cputime;
    [u_prec_modified2_coluring,flag_prec_modified2_colouring,iter_prec_modified2_colouring,absres_prec_modified2_colouring] = mygmres(A,Pmodified2_colouring,rhs,params,initialguess);
    finish_prec_modified2_colouring = cputime;

    tot_prec_modified2_colouring = tot_prec_modified2_colouring + (finish_prec_modified2_colouring - start_prec_modified2_colouring);   

end


if flag_prec_modified2_nocolouring == 0
    display(strcat('Modified Preconditioned GMRES with smoothed interpolator took:  ', num2str(finish_prec_modified2_nocolouring - start_prec_modified2_nocolouring)));
else
    display('Modified Preconditioned GMRES with smoothed interpolator did NOT converge');
end


if flag_prec_modified2_colouring == 0
    display(strcat('Modified Preconditioned GMRES with smoothed interpolator took:  ', num2str(finish_prec_modified2_colouring - start_prec_modified2_colouring)));
else
    display('Modified Preconditioned GMRES with smoothed interpolator did NOT converge');
end
