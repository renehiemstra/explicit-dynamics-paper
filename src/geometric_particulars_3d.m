    function geometry = geometric_particulars_3d(mapping, quadrature)

    % compute b-spline basis functions, up to second derivatives
    b.u = shapefuns(mapping.basis{1}, quadrature.rule{1}.points);
    b.v = shapefuns(mapping.basis{2}, quadrature.rule{2}.points);
    b.w = shapefuns(mapping.basis{3}, quadrature.rule{3}.points);

    % get coordinates
    geometry.coordinates = {tensorconstract(mapping.cpts{1}, b.u{1}, b.v{1}, b.w{1}, 2);
                            tensorconstract(mapping.cpts{2}, b.u{1}, b.v{1}, b.w{1}, 2);
                            tensorconstract(mapping.cpts{3}, b.u{1}, b.v{1}, b.w{1}, 2)};
    
    % compute gradient
    geometry.jacobian = { tensorconstract(mapping.cpts{1}, b.u{2}, b.v{1}, b.w{1}, 2), tensorconstract(mapping.cpts{1}, b.u{1}, b.v{2}, b.w{1}, 2), tensorconstract(mapping.cpts{1}, b.u{1}, b.v{1}, b.w{2}, 2);
                          tensorconstract(mapping.cpts{2}, b.u{2}, b.v{1}, b.w{1}, 2), tensorconstract(mapping.cpts{2}, b.u{1}, b.v{2}, b.w{1}, 2), tensorconstract(mapping.cpts{2}, b.u{1}, b.v{1}, b.w{2}, 2);
                          tensorconstract(mapping.cpts{3}, b.u{2}, b.v{1}, b.w{1}, 2), tensorconstract(mapping.cpts{3}, b.u{1}, b.v{2}, b.w{1}, 2), tensorconstract(mapping.cpts{3}, b.u{1}, b.v{1}, b.w{2}, 2)};

    % compute matrix of second derivatives
    for j=3:-1:1
        geometry.hessian{1,j} = tensorconstract(mapping.cpts{j}, b.u{3}, b.v{1}, b.w{1}, 2);
        geometry.hessian{2,j} = tensorconstract(mapping.cpts{j}, b.u{1}, b.v{3}, b.w{1}, 2);
        geometry.hessian{3,j} = tensorconstract(mapping.cpts{j}, b.u{1}, b.v{1}, b.w{3}, 2);
        geometry.hessian{4,j} = tensorconstract(mapping.cpts{j}, b.u{1}, b.v{2}, b.w{2}, 2);
        geometry.hessian{5,j} = tensorconstract(mapping.cpts{j}, b.u{2}, b.v{1}, b.w{2}, 2);
        geometry.hessian{6,j} = tensorconstract(mapping.cpts{j}, b.u{2}, b.v{2}, b.w{1}, 2);
    end
end