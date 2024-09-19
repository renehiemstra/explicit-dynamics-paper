% computes 'F - K * d' in a matrix-free fashion for a 2d problem
% coefs is a matrix of coefficients for 'd'
function f = assemble_standard_forcing_vector_3d(coeffs, discretization, material, forces, flag_mex)

    % initialize
    geometry = discretization.geometry;
    basis = discretization.basis;
    quadrature = discretization.quadrature;

    % compute partial derivatives of the displacment field
    displacement.ders = {tensorconstract(coeffs, basis{1}.trialfuns{2}, basis{2}.trialfuns{1}, basis{3}.trialfuns{1}, 2);
                         tensorconstract(coeffs, basis{1}.trialfuns{1}, basis{2}.trialfuns{2}, basis{3}.trialfuns{1}, 2);
                         tensorconstract(coeffs, basis{1}.trialfuns{1}, basis{2}.trialfuns{1}, basis{3}.trialfuns{2}, 2)};

    % initialize the quadrature loop
    m1 = quadrature.rule{1}.npoints;
    m2 = quadrature.rule{2}.npoints;
    m3 = quadrature.rule{3}.npoints;
    kappa = material.kappa;
    weights = {quadrature.rule{1}.weights, quadrature.rule{2}.weights, quadrature.rule{3}.weights};
    coordinates = geometry.coordinates;
    forces_eval = forces(coordinates{1}, coordinates{2}, coordinates{3});
  
    if (flag_mex==true)
        % call C function for fast evaluation of the forces
        % input is checked on matlab side
        assert(size(forces_eval, 1)==m1);
        assert(size(forces_eval, 2)==m2);
        assert(size(forces_eval, 3)==m3);
        assert(isscalar(kappa));
        % call C function
        mex_quadrature_loop_standard_3d(m1, m2, m3, ...
            kappa, ...
            quadrature.sumfact, ...
            weights, ...
            forces_eval, ...
            displacement.ders,...
            geometry.jacobian);
    else
        quadrature.sumfact = mat_quadrature_loop_standard_3d(m1, m2, m3, ...
            kappa, ...
            quadrature.sumfact, ...
            weights, ...
            forces_eval, ...
            displacement.ders,...
            geometry.jacobian);
    end

    % perform sumfactorization and yield vector
    f = tensorconstract(quadrature.sumfact{1}, basis{1}.testfuns{2}, basis{2}.testfuns{1}, basis{3}.testfuns{1}, 1) + ...
        tensorconstract(quadrature.sumfact{2}, basis{1}.testfuns{1}, basis{2}.testfuns{2}, basis{3}.testfuns{1}, 1) + ...
        tensorconstract(quadrature.sumfact{3}, basis{1}.testfuns{1}, basis{2}.testfuns{1}, basis{3}.testfuns{2}, 1) + ...
        tensorconstract(quadrature.sumfact{4}, basis{1}.testfuns{1}, basis{2}.testfuns{1}, basis{3}.testfuns{1}, 1);
    
    % invert mass matrix
    if strcmp(discretization.method, 'standard')
        f = reshape(discretization.chol \ (discretization.chol' \ f(:)), discretization.dims);
    elseif strcmp(discretization.method, 'lumped')
        f = discretization.approxinverse .* f;
    else
        error(strcat('Routine not implemented for ', discretization.method));
    end
end