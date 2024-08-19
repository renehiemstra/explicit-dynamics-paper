%% initialize
clear all; close all; clc;
addpath('../src','../mex');

% parameters
p  = 3;         % polynomial degree
ne = [2,3];     % number of elements in each direction
alpha = 0;      % parameter in [0,1] that induces nonlinearity in the mapping
k = 1;          % wavenumber in space.
lambda = 1;     % wavenumber in time.
kappa = 1;      % material diffusion coefficient
rho = 1;        % density     

flag_mex = true;

% problem specifics
problem.geometry.mapping = mapping(alpha);          % get the mapping
problem.material.kappa = kappa;                     % define the material tensor

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

%% test stiffness term of matrix-free with matrix-based

% pick random solution coefficients
n = [discretization.basis{1}.dim, discretization.basis{2}.dim];
coeffs = rand(n(1), n(2));
% stiffness term computed using matrix-free approach
f = assemble_forcing_vector(coeffs, discretization, problem.material, @(x,y) 0.0 .* x, flag_mex);

M1 = discretization.basis{1}.testfuns{1}' * ...
        (discretization.quadrature.rule{1}.weights .* ...
            discretization.basis{1}.trialfuns{1});
M2 = discretization.basis{2}.testfuns{1}' * ...
        (discretization.quadrature.rule{2}.weights .* ...
            discretization.basis{2}.trialfuns{1});
K1 = discretization.basis{1}.testfuns{2}' * ...
        (discretization.quadrature.rule{1}.weights .* ...
            discretization.basis{1}.trialfuns{2});
K2 = discretization.basis{2}.testfuns{2}' * ...
        (discretization.quadrature.rule{2}.weights .* ...
            discretization.basis{2}.trialfuns{2});

% stiffness term computed using kronecker products
f_kron = - M1 * coeffs * K2 - K1 * coeffs * M2;

% test if results are the same
assert(norm(f - f_kron) < 1e-12);

%% test with non-zero right-handside

forcing = @(x,y) sin(x).*cos(y);
x = discretization.geometry.coordinates{1};
y = discretization.geometry.coordinates{2};
wu = discretization.quadrature.rule{1}.weights;
wv = discretization.quadrature.rule{2}.weights;

% matrix-free
f = assemble_forcing_vector(coeffs, discretization, problem.material, forcing, flag_mex);

% using matrix based approach
f_kron = f_kron + (wu .* discretization.basis{1}.testfuns{1})' * (forcing(x, y) * (wv .* discretization.basis{2}.testfuns{1}));

% test if results are the same
assert(norm(f-f_kron) < 1e-12);
