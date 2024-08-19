function rule = gaussrule(basis)
    [x,w] = lgwtglobal(basis.p+1, unique(basis.kts)); m = length(x);
    rule.points = x;
    rule.weights = w;
    rule.npoints = m;
end