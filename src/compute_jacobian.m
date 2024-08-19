function A = compute_jacobian(dX, k, l)
    % matrix of first derivatives (Jacobian matrix) at quadrature point (k,l)
    A = [dX{1,1}(k, l) dX{1,2}(k, l);
         dX{2,1}(k, l) dX{2,2}(k, l)];
end