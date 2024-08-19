function funs = shapefuns(basis, u)

    % compute basis functions up to second derivatives
    uu = repmat(u', 3, 1); uu = uu(:);
    bb = spcol(basis.kts,basis.p+1,uu);

    % generate B-spline basis object
    funs = {bb(1:3:end,:), bb(2:3:end,:), bb(3:3:end,:)};
end