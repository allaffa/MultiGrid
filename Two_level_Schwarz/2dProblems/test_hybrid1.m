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

Nx = 201;
Ny = 201;
x1 = 0.0;
x2 = 2.0;
y1 = 0.0;
y2 = 2.0;

%% Problem definition

display('Construction of the problem on the fine mesh')

problem = @(x1,x2,y1,y2,Nx,Ny) fd2d_heat_steady_test05 ( x1,x2,y1,y2,Nx,Ny );

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
ndomains_x = 20;
ndomains_y = 20;

n_domains = ndomains_x * ndomains_y;

% overlapping is expressed in percentage with respect to the size of the
% single subdomain. Therefore it must be a number between 0 and 1
overlapping = 0.2;

display('Construction of the DD data structure')

[DD_data, count_contributions] = PrepareDD_Data( Nx, x1, x2, ndomains_x, Ny, y1, y2, ndomains_y, overlapping, problem );  

%% Set up multicore policy

startup

%% Construction of the Two Level Schwarz method preconditioner

P1 = @(x) Hybrid1Schwarz_Preconditioner(A, MGData, DD_data, count_contributions, x);

P2 = @(x) Hybrid2Schwarz_Preconditioner(A, MGData, DD_data, count_contributions, x);

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
maxit = 10^(4);

params(1) = tol; params(2) = maxit;  % tolerance and max iterations

initialguess = rand(size(A,1),1);

start_prec1 = cputime;
[u_prec1,flag_prec1,iter_prec1,absres_prec1]=mygmres(A, P1,rhs,params,initialguess);
finish_prec1 = cputime;

if flag_prec1 == 0
    display(strcat('Preconditioned GMRES took:  ', num2str(finish_prec1 - start_prec1)));
else
    display('Preconditioned GMRES did NOT converge');
end

start_prec2 = cputime;
[u_prec2,flag_prec2,iter_prec2,absres_prec2]=mygmres(A, P2,rhs,params,initialguess);
finish_prec2 = cputime;

if flag_prec2 == 0
    display(strcat('Preconditioned GMRES took:  ', num2str(finish_prec2 - start_prec2)));
else
    display('Preconditioned GMRES did NOT converge');
end

