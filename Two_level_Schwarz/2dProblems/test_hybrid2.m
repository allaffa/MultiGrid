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

Nx = 401;
Ny = 401;
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
ratiox = 20;
ratioy = 20;

% Relaxation technique used
relaxing = 'Jacobi';

display('Construction of the MG data structure')

% In the data structure, 1 is the index for the fine level, 2 for the
% coarse level
MGData = PrepareMGData( Nx, x1,x2, ratiox, Ny, y1, y2, ratioy, relaxing, problem );   

%% Domain Decomposition setting

% Setting the number of domains
ndomains_x = 20;
ndomains_y = 20;

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

colouring = 0;

start = cputime;
[Smoothed_Prolungation] = Hybrid2Schwarz_Smoothed_interpolation(A, MGData, DD_data, count_contributions, colouring);
finish = cputime;

Ibar_time = finish - start;

%% Construction of the Two Level Schwarz method preconditioner

P = @(x) Hybrid2Schwarz_Preconditioner(A, MGData, DD_data, count_contributions, x);

% Vassilevskj approach
Pmodified = @(x) Hybrid2Schwarz_Preconditioner_modified(A, MGData, DD_data, count_contributions, x);

colour = 0;
% Vassilevskj approach with explicitly built smoothed inteprolator
Pmodified2 = @(x) Hybrid2Schwarz_Preconditioner_modified2(A, Smoothed_Prolungation, MGData, DD_data, count_contributions, x, colour);

% Since Gmres command requires necessarily a preconditioner, I define the
% identity operator to use when I do not want to apply any preconditioner
NonP = @(x) x;

%% Definition of global variables to measure performance

global DD_time Ac_time P_time ...
    DD_time_modified Ac_time_modified P_time_modified ...
    DD_time_modified2 Ac_time_modified2 P_time_modified2 Pbar_time_modified2

DD_time = 0.0;
Ac_time = 0.0;
P_time = 0.0;

DD_time_modified = 0.0;
Ac_time_modified = 0.0;
P_time_modified = 0.0;

DD_time_modified2 = 0.0;
Ac_time_modified2 = 0.0;
P_time_modified2 = 0.0;
Pbar_time_modified2 = 0.0;


%% Call to Krylov solver 

display('Solver called');
number_runs = 100;

tot_prec = 0.0;
tot_prec_modified = 0.0;
tot_prec_modified2 = 0.0;

tol = 10^(-8);
%restart = 10; Since we are using Full GMRES actually we do not need the
%restart parameter 
maxit = 10^(3);

params(1) = tol; params(2) = maxit;  % tolerance and max iterations


for runs = 1:number_runs
    
    initialguess = rand(size(A,1),1);

    start_prec = cputime;
    [u_prec,flag_prec,iter_prec,absres_prec]=mygmres(A, P,rhs,params,initialguess);
    finish_prec = cputime;

    tot_prec = tot_prec + (finish_prec - start_prec);

    start_prec_modified = cputime;
    [u_prec_modified,flag_prec_modified,iter_prec_modified,absres_prec_modified] = mygmres(A,Pmodified,rhs,params,initialguess);
    finish_prec_modified = cputime;

    tot_prec_modified = tot_prec_modified + (finish_prec_modified - start_prec_modified);

    start_prec_modified2 = cputime;
    [u_prec_modified2,flag_prec_modified2,iter_prec_modified2,absres_prec_modified2] = mygmres(A,Pmodified2,rhs,params,initialguess);
    finish_prec_modified2 = cputime;

    tot_prec_modified2 = tot_prec_modified2 + (finish_prec_modified2 - start_prec_modified2);

end

tot_prec = tot_prec / number_runs;
tot_prec_modified = tot_prec_modified / number_runs;
tot_prec_modified2 = tot_prec_modified2 / number_runs;

DD_time = DD_time / number_runs;
Ac_time = Ac_time / number_runs;
P_time = P_time / number_runs;

DD_time_modified = DD_time_modified / number_runs;
Ac_time_modified = Ac_time_modified / number_runs;
P_time_modified = P_time_modified / number_runs;

DD_time_modified2 = DD_time_modified2 / number_runs;
Ac_time_modified2 = Ac_time_modified2 / number_runs;
P_time_modified2 = P_time_modified2 / number_runs;
Pbar_time_modified2 = Pbar_time_modified2 / number_runs;

if flag_prec == 0
    display(strcat('Preconditioned GMRES took:  ', num2str(finish_prec - start_prec)));
else
    display('Preconditioned GMRES did NOT converge');
end

if flag_prec_modified == 0
    display(strcat('Modified Preconditioned GMRES took:  ', num2str(finish_prec_modified - start_prec_modified)));
else
    display('Modified Preconditioned GMRES did NOT converge');
end

if flag_prec_modified2 == 0
    display(strcat('Modified Preconditioned GMRES with smoothed interpolator took:  ', num2str(finish_prec_modified2 - start_prec_modified2)));
else
    display('Modified Preconditioned GMRES with smoothed interpolator did NOT converge');
end