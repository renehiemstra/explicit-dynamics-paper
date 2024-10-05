function omega_max = estimate_maximum_eigenvalue_2d(discretization, problem, flag_mex)

    % dimensions of discretization
    for i=2:-1:1
        dims(i) = discretization.basis{i}.dim;
    end

    % choose solver object
    if strcmp(discretization.product, 'standard')
        % standard - L2 inner product
        solver = @(d) assemble_standard_forcing_vector_2d(d, ...
                    discretization, ...
                    problem.material, ...
                    @(x,y) 0.0 .* x, ...
                    flag_mex);
    elseif strcmp(discretization.product, 'weighted')
        % weighted L2 inner product
        solver = @(d) assemble_weighted_forcing_vector_2d(d, ...
                    discretization, ...
                    problem.material, ...
                    @(x,y) 0.0 .* x, ...
                    flag_mex);
    end

    % compute maximum eigenvalue of the stiffness matrix
    A = @(z) -reshape( solver(reshape(z, dims)), prod(dims), 1);
    omega_max = sqrt(real(eigs(A, prod(dims), 1, 'largestabs')));
end