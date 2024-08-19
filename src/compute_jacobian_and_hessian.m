function [A, DA] = compute_jacobian_and_hessian(dX, ddX, k, l)

    % matrix of first derivatives (Jacobian matrix) at quadrature point (k,l)
    A = [dX{1,1}(k, l) dX{1,2}(k, l);
         dX{2,1}(k, l) dX{2,2}(k, l)];

    % matrix of second derivatives at quadrature point (k,l)
    DA = [ddX{1,1}(k, l) ddX{1,2}(k, l);
          ddX{2,1}(k, l) ddX{2,2}(k, l);
          ddX{3,1}(k, l) ddX{3,2}(k, l)];
end