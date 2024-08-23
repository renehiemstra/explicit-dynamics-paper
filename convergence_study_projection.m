%% initialize
clear all; close all; clc;
addpath('src','mex');

% get problem specs
alpha = 0.0;        % parameter in [0,1] that induces nonlinearity in the mapping
k = 1;              % wavenumber in space.
lambda = 2;         % wavenumber in time.
kappa = 1;          % material diffusion coefficient
problem = benchmark_1_2d(k, lambda, kappa, alpha);

%% get discretization and time integrator

NE = [4, 8, 16, 32];
P = [3,4,5];
Methods = {'standard', 'standard';'standard', 'lumped'; 'weighted', 'dual'};
l2e = zeros(length(NE),length(Methods), length(P));

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
        
            % get discretization object
            [discretization, integrator] = get_discretization(problem, p, [ne,ne], 1, method, product, true);
            
            % project
            discretization.coeffs = project(discretization, @(x,y) problem.solution.displacement(x,y,0.0));

            % compute relative l2 error
            l2e(i, j, k) = evaluate_l2error(discretization, @(x,y) problem.solution.displacement(x,y,0.0), "standard");
            i = i + 1;
        end
    end
end

%% postprocessing
close all; clc;

slope = 2;
dx = 2;
dy = slope * dx;

% origin
x_a = 8;
y_a = 0.1;

x_b = x_a * dx;
y_b = y_a / dy;


% Create figure and axis
figure1 = figure;
axes1 = axes('Parent',figure1); hold(axes1,'on');
plot(NE, reshape(l2e(:,1,:), length(NE), length(P)),'Marker','o','LineWidth',2,'Parent',axes1);

plot([x_a;x_b], [y_a;y_b],'LineWidth',3, 'Color', 'b');
plot([x_a;x_b], [y_b;y_b],'LineWidth',1, 'Color', 'b');
plot([x_a;x_a], [y_a;y_b],'LineWidth',1, 'Color', 'b');

ylabel('L2-error');
xlabel('n_e');
box(axes1,'on');
hold(axes1,'off');
set(axes1,'XMinorTick','on','XScale','log','XTick', NE,...
    'YMinorTick','on','YScale','log');