%% initialize
clear all; close all; clc;
addpath('src','mex');

% parameters
p  = 3;             % polynomial degree
ne = [8,8,8];       % number of elements in each direction
Tmax = 1;           % Maximum time to simulate
flag_mex = true;   % flag to do quadrature loop in C

% pick the discretization method:
method = 'standard';    % choose: 'standard' or 'dual' or 'lumped'
product = 'standard';   % choose: 'standard' or 'weighted'

% get problem specs
alpha = 0.0;        % parameter in [0,1] that induces nonlinearity in the mapping
k = 1;              % wavenumber in space.
lambda = 2;         % wavenumber in time.
kappa = 1;          % material diffusion coefficient
problem = benchmark_1_3d(k, lambda, kappa, alpha);

%% get discretization and time integrator

% get discretization object
[discretization, integrator] = get_discretization_3d(problem, p, ne, Tmax, method, product, flag_mex);

% Runge-Kutta time integration
discretization = solve_3d(problem, discretization, integrator, flag_mex);

%% postprocessing

% evaluate the mapping
plot_2d_contours(problem, discretization, [20,20], 0.5)

% compute relative l2 error
l2e = evaluate_l2error_3d(discretization, @(x,y,z) problem.solution.displacement(x,y,z,Tmax), "standard");
fprintf('\n\nThe relative L2 error is: %0.2e\n', l2e);