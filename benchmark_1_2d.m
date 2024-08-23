function problem = benchmark_1_2d(k, lambda, kappa, alpha)

    % problem specifics
    problem.geometry.mapping = mapping(alpha);  % get the mapping
    problem.material.kappa = kappa;             % define the material tensor

    % exact benchmark solution
    problem.solution.displacement = @(x,y,t) sin(k*pi*x) .* sin(k*pi*y) .* cos(pi*lambda*t);  
    problem.solution.velocity = @(x,y,t) -lambda*pi * sin(pi*k*x) .* sin(pi*k*y) .* sin(pi*lambda*t); % velocity
    problem.solution.forcing  = @(x,y,t) pi^2 * (2 * kappa * k^2 - lambda^2) * problem.solution.displacement(x,y,t);

end