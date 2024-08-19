%% initialize
clear all; close all; clc;
addpath('src','mex');

% parameters
p  = 3;             % polynomial degree
ne = [128,128];     % number of elements in each direction
alpha = 0.0;        % parameter in [0,1] that induces nonlinearity in the mapping
k = 2;              % wavenumber in space.
lambda = 2;         % wavenumber in time.
kappa = 1;          % material diffusion coefficient
rho = 1;            % density
Tmax = 1.0;         % Maximum time to simulate
flag_mex = true;    % flag to do quadrature loop in C

% pick the discretization method:
discretization.method = 'dual'; % choose: 'standard' or 'dual' or 'lumped'
discretization.product = 'weighted'; % choose: 'standard' or 'weighted'

% problem specifics
problem.geometry.mapping = mapping(alpha);  % get the mapping
problem.material.kappa = kappa;             % define the material tensor
% exact benchmark solution
problem.solution.displacement = @(x,y,t) sin(k*pi*x) .* sin(k*pi*y) .* cos(pi*lambda*t);  
problem.solution.velocity = @(x,y,t) -lambda*pi * sin(pi*k*x) .* sin(pi*k*y) .* sin(pi*lambda*t); % velocity
problem.solution.forcing  = @(x,y,t) - pi^2 * (lambda^2 - 2 * kappa * k^2) * problem.solution.displacement(x,y,t);

% construct basis u-dir
discretization.basis{1}.p = p;
discretization.basis{1}.kts = knotvector(discretization.basis{1}.p, 0.0 ,1.0 ,ne(1));
discretization.basis{1}.eval.trialfuns = @(u) trialfuns(discretization.basis{1}, u);
discretization.basis{1}.eval.testfuns = @(u) testfuns(discretization.basis{1}, u);

% construct basis v-dir
discretization.basis{2}.p = p;
discretization.basis{2}.kts = knotvector(discretization.basis{2}.p, 0.0 ,1.0 ,ne(2));
discretization.basis{2}.eval.trialfuns = @(v) trialfuns(discretization.basis{2}, v);
discretization.basis{2}.eval.testfuns = @(v) testfuns(discretization.basis{2}, v);

% generate quadrature object
discretization.quadrature.rule{1} = gaussrule(discretization.basis{1});
discretization.quadrature.rule{2} = gaussrule(discretization.basis{2});
m = [discretization.quadrature.rule{1}.npoints, discretization.quadrature.rule{2}.npoints];
% allocate space for performing sum factorization
discretization.quadrature.sumfact = {zeros(m(1), m(2)), zeros(m(1), m(2)), zeros(m(1), m(2))};

% precompute geometric particulars, including first and second derivatives
discretization.geometry = geometric_particulars(problem.geometry.mapping, discretization.quadrature);

% pre-compute b-spline basis functions
% basis in u-dir
discretization.basis{1}.trialfuns = discretization.basis{1}.eval.trialfuns(discretization.quadrature.rule{1}.points);
discretization.basis{1}.testfuns = discretization.basis{1}.eval.testfuns(discretization.quadrature.rule{1}.points);
discretization.basis{1}.dim = size(discretization.basis{1}.testfuns{1},2);
% basis in v-dir
discretization.basis{2}.trialfuns = discretization.basis{2}.eval.trialfuns(discretization.quadrature.rule{2}.points);
discretization.basis{2}.testfuns = discretization.basis{2}.eval.testfuns(discretization.quadrature.rule{2}.points);
discretization.basis{2}.dim = size(discretization.basis{2}.testfuns{1},2);

% compute the mass matrix
discretization.massmatrix{1} = mass_matrix(discretization.basis{1}, discretization.quadrature.rule{1});
discretization.massmatrix{2} = mass_matrix(discretization.basis{2}, discretization.quadrature.rule{2});

if strcmp(discretization.product, 'standard')
    if strcmp(discretization.method, 'standard')
        discretization.mass = assemble_mass_matrix(discretization, flag_mex);
        discretization.chol = chol(discretization.mass);
    elseif strcmp(discretization.method, 'lumped')
        discretization.approxinverse = reshape(1 ./ sum(assemble_mass_matrix(discretization, flag_mex), 2), discretization.basis{1}.dim, discretization.basis{2}.dim);
    else
        error(strcat('Not implemented for ', discretization.method));
    end
elseif strcmp(discretization.product, 'weighted')
    if strcmp(discretization.method, 'standard')
        discretization.chol{1} = chol(discretization.massmatrix{1});
        discretization.chol{2} = chol(discretization.massmatrix{2});
    elseif strcmp(discretization.method, 'lumped')
        discretization.approxinverse{1} = 1 ./ sum(discretization.massmatrix{1}, 2);
        discretization.approxinverse{2} = 1 ./ sum(discretization.massmatrix{2}, 2);
    elseif strcmp(discretization.method, 'dual')
        discretization.approxinverse{1} = approximate_l2_inverse(discretization.basis{1}.p, discretization.basis{1}.kts, true, true);
        discretization.approxinverse{2} = approximate_l2_inverse(discretization.basis{2}.p, discretization.basis{2}.kts, true, true);
    else
        error(strcat('Not implemented for ', discretization.method));
    end
else 
    error('Choose <standard> or <weighted> as discretization.product.');
end

%% Setup time integrator

% basic Runge-Kutta integrator setup
integrator = setup_rk_method(discretization);

% final time and safety factor applied to critical time-step
integrator.tmax = Tmax;
integrator.safety_factor = 0.3;

% compute the maximum eigenfrequency
omega_max = estimate_maximum_eigenvalue(discretization, problem);

% compute the critical time-step
integrator.dtcrit = integrator.cmax / omega_max;

% compute time-step
nT = ceil(integrator.tmax / (integrator.safety_factor * integrator.dtcrit));
integrator.dt = integrator.tmax / nT;

%% Runge-Kutta time integration

% get initial quess using L2-projection
d = project(discretization, @(x,y) problem.solution.displacement(x,y,0.0));
v = project(discretization, @(x,y) problem.solution.velocity(x,y,0.0));

% choose solver object
if strcmp(discretization.product, 'standard')
    % standard - L2 inner product
    solver = @(d, v, t) assemble_standard_forcing_vector(d, discretization, problem.material, @(x, y) problem.solution.forcing(x, y, t), flag_mex);
elseif strcmp(discretization.product, 'weighted')
    % weighted L2 inner product
    solver = @(d, v, t) assemble_weighted_forcing_vector(d, discretization, problem.material, @(x, y) problem.solution.forcing(x, y, t), flag_mex);
end

T = integrator.tmax;
dt = integrator.dt;
i = 1;
for t=dt:dt:T
    
    % status
    fprintf('time integrator, p = %d, mesh = %d x %d, time-step %d%s%d \n', p, ne(1), ne(2),  i, '/', integrator.tmax/integrator.dt);
    
    % compute new displacement, velocity and acceleration
    [d, v] = integrator.solver(d, v, t, dt, solver);

    % increase counter
    i = i+1;
end


%% postprocessing

% d = project(discretization, @(x,y) problem.solution.displacement(x,y,0.0));
% v = project(discretization, @(x,y) problem.solution.velocity(x,y,0.0));

discretization.coeffs = d;

% evaluate the mapping
gran = [20,20];
[x, y] = evaluate_mapping(problem.geometry.mapping, gran);
u_bar  = evaluate_field(discretization, gran);

% compute relative l2 error
evaluate_l2error(discretization, @(x,y) problem.solution.displacement(x,y,0.0), "relative")

% plot mesh
surf(x, y, u_bar);