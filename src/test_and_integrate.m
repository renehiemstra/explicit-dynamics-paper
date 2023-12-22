function F = test_and_integrate(fun, p, kts_r, kts_h, C)

% compute L^2 error in the mode shapes
[x_r,w_r] = lgwtglobal(p+1,unique(kts_r));
[x_h,w_h] = lgwtglobal(p+1,unique(kts_h));

% compute basis functions at quadrature points
N_r = spcol(kts_r,p+1,x_r);
N_h = spcol(kts_h,p+1,x_h);

f = (N_r' * (fun(x_r, x_h') .* ((x_r .* w_r) * w_h')) * N_h)';
F = C' * f(:);
end