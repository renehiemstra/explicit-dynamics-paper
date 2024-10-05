%% initialize
clear all; close all; clc;
addpath('src','mex');

% parameters
Tmax = 1.0;         % Maximum time to simulate
flag_mex = true;    % flag to do quadrature loop in C

% get problem specs
alpha = 0.0;        % parameter in [0,1] that induces nonlinearity in the mappmethoding
kwave = 4;          % wavenumber in space.
lambda = 4;         % wavenumber in time.
kappa = 1;          % material diffusion coefficient
problem = benchmark_1_2d(kwave, lambda, kappa, alpha);

%% get discretization and time integrator

NE = [16,32,64,128];
P = [3,4,5];
Methods = {'standard', 'standard';'standard', 'lumped'; 'weighted', 'dual'};
l2e = zeros(length(NE),length(Methods), length(P)); time = zeros(length(NE),length(Methods), length(P));

% loop over discretization approaches
for k=1:length(P)
    % polynomial degree
    p = P(k);
    
    for j=1:size(Methods, 1)
        
        product = Methods{j, 1};
        method  = Methods{j, 2};
    
        % loop over elements
        i = 1;
        for ne = NE
        
            % print case
            fprintf('case: %s-%s, p = %d, n_e = %d\n', method, product, p, ne);
        
            % measure runtime
            tic;
                % get discretization object
                [discretization, integrator] = get_discretization_2d(problem, p, [ne,ne], Tmax, method, product, flag_mex);
                
                % Runge-Kutta time integration
                discretization = solve_2d(problem, discretization, integrator, flag_mex);
            time(i, j, k) = toc;
        
            % compute relative l2 error
            l2e(i, j, k) = evaluate_l2error_2d(discretization, @(x,y) problem.solution.displacement(x,y,Tmax), "relative");
            i = i + 1;
        end
    end
end

filename = sprintf('study-2d-alpha%d-lambda%d-k%d', 10*alpha, lambda, kwave);
save(filename, "NE","P","time","l2e");