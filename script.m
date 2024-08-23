%% initialize
clear all; close all; clc;
addpath('src','mex');

% parameters
p  = 5;             % polynomial degree
ne = [64,64];       % number of elements in each direction
Tmax = 1;         % Maximum time to simulate
flag_mex = true;    % flag to do quadrature loop in C

% pick the discretization method:
method = 'standard'; % choose: 'standard' or 'dual' or 'lumped'
product = 'standard'; % choose: 'standard' or 'weighted'

% get problem specs
alpha = 0.0;        % parameter in [0,1] that induces nonlinearity in the mapping
k = 1;              % wavenumber in space.
lambda = 2;         % wavenumber in time.
kappa = 1;          % material diffusion coefficient
problem = benchmark_1_2d(k, lambda, kappa, alpha);

%% get discretization and time integrator

% get discretization object
[discretization, integrator] = get_discretization(problem, p, ne, Tmax, method, product, flag_mex);

MM = system_matrix(discretization.basis{1}, discretization.quadrature.rule{1}, 1, 1);
KK = system_matrix(discretization.basis{1}, discretization.quadrature.rule{1}, 2, 2);

M = kron(MM, MM);
K = -kron(MM, KK) - kron(KK, MM);


m = discretization.basis{1}.dim;
n = discretization.basis{2}.dim;
discretization.coeffs = rand(m , n);
Z = M \ (K * discretization.coeffs(:));

f = assemble_standard_forcing_vector(discretization.coeffs, discretization, problem.material, @(x, y) 0.0*x, flag_mex);




% Runge-Kutta time integration
discretization = solve(problem, discretization, integrator, flag_mex);

%% postprocessing

% evaluate the mapping
plot_2d_contours(problem, discretization, [20,20])

% compute relative l2 error
l2e = evaluate_l2error(discretization, @(x,y) problem.solution.displacement(x,y,Tmax), "standard");
fprintf('\n\nThe relative L2 error is: %0.2e\n', l2e);
