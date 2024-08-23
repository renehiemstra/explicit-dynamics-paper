%% initialize
clear all; close all; clc;
addpath('src','mex');

% parameters
Tmax = 1.0;         % Maximum time to simulate
flag_mex = true;    % flag to do quadrature loop in C

% get problem specs
alpha = 0.0;        % parameter in [0,1] that induces nonlinearity in the mappmethoding
k = 2;              % wavenumber in space.
lambda = 2;         % wavenumber in time.
kappa = 1;          % material diffusion coefficient
problem = benchmark_1_2d(k, lambda, kappa, alpha);

%% get discretization and time integrator

NE = [16, 32, 64, 128];
P = [5];
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
                [discretization, integrator] = get_discretization(problem, p, [ne,ne], Tmax, method, product, flag_mex);
                
                % Runge-Kutta time integration
                discretization = solve(problem, discretization, integrator, flag_mex);
            time(i, j, k) = toc;
        
            % compute relative l2 error
            l2e(i, j, k) = evaluate_l2error(discretization, @(x,y) problem.solution.displacement(x,y,Tmax), "relative");
            i = i + 1;
        end
    end
end

%% postprocessing
close all; clc;

% Create figure and axis
figure1 = figure;

axes1 = axes('Parent',figure1); hold(axes1,'on');
plot(NE, l2e,'Marker','o','LineWidth',2,'Parent',axes1);
createtriangle(75, 10^-7.4, 4, 0.2, 'b');
createtriangle(75, 10^-2.3, 2, 0.2, 'r');

ylabel('L2-error');
xlabel('n_e');
box(axes1,'on');
hold(axes1,'off');
set(axes1,'XMinorTick','on','XScale','log','XTick', NE,...
    'YMinorTick','on','YScale','log');


figure; loglog(NE, time, '-o', 'LineWidth', 2);
xlabel('n_e'); ylabel('time'); xticks(NE);
figure; loglog(time, l2e, '-o', 'LineWidth', 2);
ylabel('L2-error'); xlabel('time');