function funs = testfuns(basis, u)

    % Extraction operator:
    %   incorporates homogeneous boundary conditions
    E = extraction_operator(constraint_matrix(basis.p, basis.kts, 1:1, 1:1));

    % compute basis functions up to second derivatives
    uu = repmat(u', 3, 1); uu = uu(:);
    bb = spcol(basis.kts,basis.p+1,uu);

    % generate B-spline basis object
    funs = {bb(1:3:end,:) * E, bb(2:3:end,:) * E, bb(3:3:end,:) * E};
end