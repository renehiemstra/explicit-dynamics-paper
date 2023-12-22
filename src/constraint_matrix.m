function A = constraint_matrix(p, kts, v1, v2)
    c_left  = spcol(kts, p+1, kts(1) * ones(p+1,1));    % left boundary
    c_right = spcol(kts, p+1, kts(end) * ones(p+1,1));  % right boundary
    A = [c_left(v1,:); c_right(v2,:)];
end