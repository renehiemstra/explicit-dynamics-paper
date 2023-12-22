function l2e_rel = l2error_circular_membrane(p, kts_r,kts_h,u,u_h)

% quadrature rules
[x_r,w_r] = lgwtglobal(2*p+1,unique(kts_r));
[x_h,w_h] = lgwtglobal(2*p+1,unique(kts_h(p+1:end-p)));
W_r = diag(w_r.*x_r);
W_h = diag(w_h);

% error at the quadrature points
U = u(x_r,x_h');
e = U - u_h(x_r,x_h');

% L2 error
l2u = sqrt(sum(sum(W_r * (U.^2) * W_h)));
l2e = sqrt(sum(sum(W_r * (e.^2) * W_h)));
l2e_rel = l2e / l2u; 
end