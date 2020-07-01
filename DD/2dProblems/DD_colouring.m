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


ndomains_x = 20;
ndomains_y = 20;

n_domains = ndomains_x * ndomains_y;

overlapping = 0.0;

problem = @(x1,x2,y1,y2,Nx,Ny) fd2d_heat_steady_test02 ( x1,x2,y1,y2,Nx,Ny );

[DD_data, count_contributions] = PrepareDD_Data( Nx, x1, x2, ndomains_x, Ny, y1, y2, ndomains_y, overlapping, problem );  

[W, DD_data] = colouring(DD_data, Nx*Ny, ndomains_x);
