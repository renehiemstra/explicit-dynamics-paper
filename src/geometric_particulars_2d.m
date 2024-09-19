    function geometry = geometric_particulars_2d(mapping, quadrature)

    % compute b-spline basis functions, up to second derivatives
    b.u = shapefuns(mapping.basis{1}, quadrature.rule{1}.points);
    b.v = shapefuns(mapping.basis{2}, quadrature.rule{2}.points);

    % get coordinates
    geometry.coordinates = {tensorconstract(mapping.cpts{1}, b.u{1}, b.v{1}, 2);
                            tensorconstract(mapping.cpts{2}, b.u{1}, b.v{1}, 2)};
    
    % compute gradient
    geometry.jacobian = { tensorconstract(mapping.cpts{1}, b.u{2}, b.v{1}, 2), tensorconstract(mapping.cpts{1}, b.u{1}, b.v{2}, 2);
                          tensorconstract(mapping.cpts{2}, b.u{2}, b.v{1}, 2), tensorconstract(mapping.cpts{2}, b.u{1}, b.v{2}, 2)};

    % compute matrix of second derivatives
    geometry.hessian = { tensorconstract(mapping.cpts{1}, b.u{3}, b.v{1}, 2), tensorconstract(mapping.cpts{2}, b.u{3}, b.v{1}, 2);
                         tensorconstract(mapping.cpts{1}, b.u{1}, b.v{3}, 2), tensorconstract(mapping.cpts{2}, b.u{1}, b.v{3}, 2);
                         tensorconstract(mapping.cpts{1}, b.u{2}, b.v{2}, 2), tensorconstract(mapping.cpts{2}, b.u{2}, b.v{2}, 2)};
end