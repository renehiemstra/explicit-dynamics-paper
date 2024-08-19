function geometry = geometric_particulars(mapping, quadrature)

    % compute b-spline basis functions, up to second derivatives
    b.u = shapefuns(mapping.basis{1}, quadrature.rule{1}.points);
    b.v = shapefuns(mapping.basis{2}, quadrature.rule{2}.points);

    % get coordinates
    geometry.coordinates = {b.u{1} * mapping.cpts{1} * b.v{1}', b.u{1} * mapping.cpts{2} * b.v{1}'};

    % compute gradient
    geometry.jacobian = { b.u{2} * mapping.cpts{1} * b.v{1}', b.u{1} * mapping.cpts{1} * b.v{2}';
                          b.u{2} * mapping.cpts{2} * b.v{1}', b.u{1} * mapping.cpts{2} * b.v{2}'};

    % compute matrix of second derivatives
    geometry.hessian = {b.u{3} * mapping.cpts{1} * b.v{1}', b.u{3} * mapping.cpts{2} * b.v{1}';
                        b.u{1} * mapping.cpts{1} * b.v{3}', b.u{1} * mapping.cpts{2} * b.v{3}';
                        b.u{2} * mapping.cpts{1} * b.v{2}', b.u{2} * mapping.cpts{2} * b.v{2}'};
end