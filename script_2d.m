%% initialize
clear all; close all; clc;
addpath('src','mex','altmany-export_fig-v3/');

% parameters
p  = 3;             % polynomial degree
ne = [64, 64];      % number of elements in each direction
Tmax = 1;           % Maximum time to simulate
flag_mex = true;    % flag to do quadrature loop in C

% pick the discretization method:
method = 'dual'; % choose: 'standard' or 'dual' or 'lumped'
product = 'weighted'; % choose: 'standard' or 'weighted'

% get problem specs
alpha = 0;        % parameter in [0,1] that induces nonlinearity in the mapping
k = 4;              % wavenumber in space.
lambda = 4;         % wavenumber in time.
kappa = 1;          % material diffusion coefficient
problem = benchmark_1_2d(k, lambda, kappa, alpha);

%% get discretization and time integrator

% get discretization object
[discretization, integrator] = get_discretization_2d(problem, p, ne, Tmax, method, product, flag_mex);

% Runge-Kutta time integration
discretization = solve_2d(problem, discretization, integrator, flag_mex);

%% postprocessing

% evaluate the mapping
plot_2d_contours(problem, discretization, [100,100]);

% compute relative l2 error
l2e = evaluate_l2error_2d(discretization, @(x,y) problem.solution.displacement(x,y,Tmax), "standard");
fprintf('\n\nThe relative L2 error is: %0.2e\n', l2e);
