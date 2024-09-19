function problem = benchmark_1_2d(k, lambda, kappa, alpha)

    % problem specifics
    problem.geometry.mapping = mapping(alpha);  % get the mapping
    problem.material.kappa = kappa;             % define the material tensor

    % exact benchmark solution
    problem.solution.displacement = @(x,y,t) sin(k*pi*x) .* sin(k*pi*y) .* cos(pi*lambda*t);  
    problem.solution.velocity = @(x,y,t) -lambda*pi * sin(pi*k*x) .* sin(pi*k*y) .* sin(pi*lambda*t); % velocity
    problem.solution.forcing  = @(x,y,t) pi^2 * (2 * kappa * k^2 - lambda^2) * problem.solution.displacement(x,y,t);

end

function map = mapping(alpha)
% generate a B-spline mapping of the unit-square
% with a non-linear mapping parameterized via 'alpha'

    % spline discretization
    % both directions the same discretization is used
    p = 3;
    kts = knotvector(p, 0.0 ,1.0 ,1);

    % contol points of the mapping
    gp = grevillepts(p, kts);
    [ycpts, xcpts] = meshgrid(gp,gp); 
    ycpts(2,2:3) = ycpts(2,2:3) + alpha;
    ycpts(3,2:3) = ycpts(3,2:3) - alpha;
    xcpts(2:3,2) = xcpts(2:3,2) + alpha;
    xcpts(2:3,3) = xcpts(2:3,3) - alpha;

    % B-spline mapping
    % x-dir
    map.basis{1}.p = p;
    map.basis{1}.kts = kts;
    map.cpts{1} = xcpts;
    % y-dir
    map.basis{2}.p = p;
    map.basis{2}.kts = kts;
    map.cpts{2} = ycpts;
end