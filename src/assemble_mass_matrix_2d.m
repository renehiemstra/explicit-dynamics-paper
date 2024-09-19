% computes 'F - K * d' in a matrix-free fashion for a 2d problem
% coefs is a matrix of coefficients for 'd'
function M = assemble_mass_matrix_2d(discretization, flag_mex)

    % initialize
    geometry = discretization.geometry;
    basis = discretization.basis;
    quadrature = discretization.quadrature;

    % initialize the quadrature loop
    m1 = quadrature.rule{1}.npoints;
    m2 = quadrature.rule{2}.npoints;
    weights = {quadrature.rule{1}.weights, quadrature.rule{2}.weights};
  
    if (flag_mex==true)
        % call C function for fast evaluation of the mass
        % input is checked on matlab side
        assert(isscalar(m1) && round(m1)==m1);
        assert(isscalar(m2) && round(m2)==m2);
        % call C function
        mex_quadrature_loop_standard_mass_2d(m1, m2, ...
            quadrature.sumfact, ...
            weights, ...
            geometry.jacobian);
    else
        quadrature.sumfact = mat_quadrature_loop_standard_mass_2d(m1, m2, ...
            quadrature.sumfact, ...
            weights, ...
            geometry.jacobian);
    end

    % compute tensor product test and trial functions
    testfuns  = kron(sparse(basis{2}.testfuns{1}), sparse(basis{1}.testfuns{1}));
    trialfuns = kron(sparse(basis{2}.trialfuns{1}), sparse(basis{1}.trialfuns{1}));

    % compute mass matrix
    M = testfuns' * (quadrature.sumfact{1}(:) .* trialfuns);

end