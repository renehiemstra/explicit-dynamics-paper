function coeffs = project_3d(discretization, fun)

    % initialize
    coordinates = discretization.geometry.coordinates;
    for i=3:-1:1
        x{i} = coordinates{i};
        w{i} = discretization.quadrature.rule{i}.weights;
        b{i} = discretization.basis{i}.testfuns{1};
        m{i} = discretization.massmatrix{i};
    end

    % compute inner product
    l2product.fun_and_testfuns = tensorconstract(fun(x{1}, x{2}, x{3}), w{1} .* b{1}, w{2} .* b{2}, w{3} .* b{3}, 1);
    
    % perform projection onto b-spline basis
    coeffs = tensorconstract(l2product.fun_and_testfuns, inv(m{1}), inv(m{2}), inv(m{3}), 1);
end