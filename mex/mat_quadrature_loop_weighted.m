% perform the quadrature loop
function sumfact = mat_quadrature_loop_weighted(m1, m2, kappa, ...
        sumfact, ...
        weights, ...
        forcing, ...
        displacement_ders,...
        jacobian,...
        hessian)

    % for each quadrature point
    for l=1:m2
        for k=1:m1
            % compute jacobian matrix of first and second derivatives
            % @ quadrature point (k,l)
            [A,DA] = compute_jacobian_and_hessian(jacobian, hessian, k, l);

            % compute the determinant of the jacobian matrix
            % @ quadrature point (k,l)
            c = det(A);

            % retreive quadrature weight
            w = weights{1}(k) * weights{2}(l);

            % compute the gradient of the determinant of the jacobian matrix
            % @ quadrature point (k,l) A' \ 
            c_grad = A' \ [DA(1,1) * A(2,2) + A(1,1) * DA(3,2) - DA(3,1) * A(2,1) - A(1,2) * DA(1,2);
                           DA(3,1) * A(2,2) + A(1,1) * DA(2,2) - DA(2,1) * A(2,1) - A(1,2) * DA(3,2)];
            
            % compute the gradient of the displacement field
            % @ quadrature point (k,l)
            d_grad = A' \ [displacement_ders{1}(k,l); displacement_ders{2}(k,l)];
            
            % heat flux @ quadrature point (k,l)
            q = kappa * d_grad;
            
            % first two elements of 'g' relate to the regular gradient.
            % the third element of 'g' relates to taking the gradient of
            % the (1/c) term and the addition of the forcing term.
            g = [-A \ q; 
                (dot(q, c_grad) / c) + forcing(k,l)];
            
            % add results to sumfactorization objects.
            % the terms in g are weighted by the quaderature weight.
            sumfact{1}(k,l) = g(1) * w; 
            sumfact{2}(k,l) = g(2) * w;
            sumfact{3}(k,l) = g(3) * w;
        end
    end
end