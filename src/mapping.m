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