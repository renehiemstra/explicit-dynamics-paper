function coeffs = project_2d(discretization, fun)

    % initialize
    coordinates = discretization.geometry.coordinates;
    x = coordinates{1}; y = coordinates{2};
    wu = discretization.quadrature.rule{1}.weights; wv = discretization.quadrature.rule{2}.weights;
    bu = discretization.basis{1}.testfuns{1}; bv = discretization.basis{2}.testfuns{1};

    % compute inner product
    l2product.fun_and_testfuns = tensorconstract(fun(x, y), wu .* bu, wv .* bv, 1);
    
    % perform projection onto b-spline basis
    coeffs = discretization.massmatrix{1} \ l2product.fun_and_testfuns / discretization.massmatrix{2};
end