function omega_max = estimate_maximum_eigenvalue(discretization, problem)

    % dimensions of discretization
    dims = [discretization.basis{1}.dim, discretization.basis{2}.dim];

    % choose solver object
    if strcmp(discretization.product, 'standard')
        % standard - L2 inner product
        solver = @(d) assemble_standard_forcing_vector(d, ...
                    discretization, ...
                    problem.material, ...
                    @(x,y) 0.0 .* x, ...
                    true);
    elseif strcmp(discretization.product, 'weighted')
        % weighted L2 inner product
        solver = @(d) assemble_weighted_forcing_vector(d, ...
                    discretization, ...
                    problem.material, ...
                    @(x,y) 0.0 .* x, ...
                    true);
    end


    % compute maximum eigenvalue of the stiffness matrix
    A = @(z) -reshape( solver(reshape(z, dims(1), dims(2))), prod(dims), 1);
    omega_max = sqrt(eigs(A, prod(dims), 1, 'largestabs'));
end