%% initialize
clear; close all; clc;
addpath('../src','../mex');

m1 = 5; m2 = 10;
kappa = 1.0;
sumfact_mex = {zeros(m1,m2), zeros(m1,m2), zeros(m1,m2)};
sumfact_mat = {zeros(m1,m2), zeros(m1,m2), zeros(m1,m2)};
weights = {rand(m1,1), rand(m2,1)};
forcing = rand(m1,m2);
displacement_ders = {rand(m1,m2), rand(m1,m2)};
jacobian = {rand(m1,m2)+ones(m1,m2), rand(m1,m2); rand(m1,m2), rand(m1,m2)+ones(m1,m2)};
hessian  = {rand(m1,m2), rand(m1,m2); rand(m1,m2), rand(m1,m2); rand(m1,m2), rand(m1,m2)};


% call C routine
mex_start = tic;
mex_quadrature_loop(m1, m2, kappa, ...
    sumfact_mex, ...
    weights, ...
    forcing, ...
    displacement_ders,...
    jacobian,...
    hessian);
mex_time = toc(mex_start);

% call matlab routine
mat_start = tic;
sumfact_mat = quadrature_loop(m1, m2, kappa, ...
        sumfact_mat, ...
        weights, ...
        forcing, ...
        displacement_ders,...
        jacobian,...
        hessian);
mat_time = toc(mat_start);

fprintf("Speedup: %f\n", mat_time / mex_time);

% check that the results from C and Matlab are the same
assert(norm(sumfact_mex{1}-sumfact_mat{1})/(m1*m2) < 1e-8);
assert(norm(sumfact_mex{2}-sumfact_mat{2})/(m1*m2) < 1e-8);
assert(norm(sumfact_mex{3}-sumfact_mat{3})/(m1*m2) < 1e-8);
