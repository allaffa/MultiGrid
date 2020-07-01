close all
clear
clc

%% Inclusion of the paths to use the tools from GMG and DD

addpath('../../StartUpSetting')
addpath('../../1dProblems')
addpath('../../1dProblems/tools')
addpath('../../GMG/1dProblems/tools')
addpath('../../DD/1dProblems/tools')
addpath('../../LinearSolvers')
addpath('../Preconditioners')

%% For the construction of the MG and DD data structures we do't need to measure the performance, so I use all the cores possible for this

maxNumCompThreads('automatic')

%% Domain setting

N = 301;
x1 = 0.0;
x2 = 1.0;

%% Problem definition

problem = @(N,x1,x2 ) fd1d_bvp_test03 ( N,x1,x2 );

%Computation of exact solution
[A,rhs] = problem ( N,x1,x2 );
exact = A\rhs;

%% Multigrid setting

% level of coarsening for the coarse mesh
Hratio = 10;

% Relaxation technique used
relaxing = 'Jacobi';

% In the data structure, 1 is the index for the fine level, 2 for the
% coarse level
MGData = PrepareMGData( N, x1,x2, Hratio, relaxing, problem );   

%% Domain Decomposition setting

% Setting the number of domains
n_domains = 50;

% overlapping is expressed in percentage with respect to the size of the
% single subdomain. Therefore it must be a number between 0 and 1
overlapping = 0.2;

[DD_data, count_contributions] = PrepareDD_Data( N, x1,x2, n_domains, overlapping, problem );   

%% Set up multicore policy

startup

%% Explicit construction of the smoothed interpolation

display( 'Construction of the smoothed preconditioner' )

colouring = 0;

[Smoothed_Prolungation] = Hybrid2Schwarz_Smoothed_interpolation(A, MGData, DD_data, count_contributions, colouring);

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

tol = 10^(-8);
%restart = 10; Since we are using Full GMRES actually we do not need the
%restart parameter 
maxit = 10^(3);

params(1) = tol; params(2) = maxit;  % tolerance and max iterations

initialguess = rand(size(A,1),1);

start_prec = cputime;
[u_prec,flag_prec,iter_prec,absres_prec]=mygmres(A, P,rhs,params,initialguess);
finish_prec = cputime;

if flag_prec == 0
    display(strcat('Preconditioned GMRES took:  ', num2str(finish_prec - start_prec)));
else
    display('Preconditioned GMRES did NOT converge');
end

start_prec_modified = cputime;
[u_prec_modified,flag_prec_modified,iter_prec_modified,absres_prec_modified] = mygmres(A,Pmodified,rhs,params,initialguess);
finish_prec_modified = cputime;

if flag_prec_modified == 0
    display(strcat('Modified Preconditioned GMRES took:  ', num2str(finish_prec_modified - start_prec_modified)));
else
    display('Modified Preconditioned GMRES did NOT converge');
end

start_prec_modified2 = cputime;
[u_prec_modified2,flag_prec_modified2,iter_prec_modified2,absres_prec_modified2] = mygmres(A,Pmodified2,rhs,params,initialguess);
finish_prec_modified2 = cputime;

if flag_prec_modified2 == 0
    display(strcat('Modified Preconditioned GMRES with smoothed iterpolator took:  ', num2str(finish_prec_modified2 - start_prec_modified2)));
else
    display('Modified Preconditioned GMRES with smoothed iterpolator did NOT converge');
end
 
start_noprec = cputime;
[u_noprec,flag_noprec,iter_noprec,absres_noprec] = mygmres(A,NonP,rhs,params,initialguess);
finish_noprec = cputime;
 
if flag_noprec == 0
    display(strcat('Non-preconditioned GMRES took:  ', num2str(finish_noprec - start_noprec)));
else
     display('Non Preconditioned GMRES did NOT converge');
end