function l2e = evaluate_l2error_2d(discretization, solution, flag)

    % initialize
    coordinates = discretization.geometry.coordinates;
    jacobian = discretization.geometry.jacobian;
    x = coordinates{1}; y = coordinates{2};
    wu = discretization.quadrature.rule{1}.weights; wv = discretization.quadrature.rule{2}.weights;

    % compute error at array of quadrature points
    vol = compute_det(jacobian);
    s = solution(x, y);
    e = s - tensorconstract(discretization.coeffs, ...
                discretization.basis{1}.trialfuns{1}, ...
                discretization.basis{2}.trialfuns{1}, ...
                2);
    
    % compute l2-error by contracting with quadrature weight vectors
    l2e = sqrt(tensorconstract((e.^2) .* vol, wu, wv, 1));

    % compute relative error
    if strcmp(flag, "relative")
        l2e = l2e / sqrt(tensorconstract((s.^2) .* vol, wu, wv, 1));
    end
end