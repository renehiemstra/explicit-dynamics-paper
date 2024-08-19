% perform the quadrature loop
function sumfact = mat_quadrature_loop_standard(m1, m2, kappa, ...
        sumfact, ...
        weights, ...
        forcing, ...
        displacement_ders,...
        jacobian)

    % for each quadrature point
    for l=1:m2
        for k=1:m1
            % compute jacobian matrix of first and second derivatives
            % @ quadrature point (k,l)
            A = compute_jacobian(jacobian, k, l);

            % compute the determinant of the jacobian matrix
            % @ quadrature point (k,l)
            c = det(A);

            % retreive quadrature weight
            w = weights{1}(k) * weights{2}(l) * c;

            % compute the gradient of the displacement field
            % @ quadrature point (k,l)
            d_grad = A' \ [displacement_ders{1}(k,l); displacement_ders{2}(k,l)];
            
            % heat flux @ quadrature point (k,l)
            q = kappa * d_grad;
            
            % first two elements of 'g' relate to the regular gradient.
            % the third element of 'g' relates to taking the gradient of
            % the (1/c) term and the addition of the forcing term.
            g = [-A \ q; 
                forcing(k,l)];
            
            % add results to sumfactorization objects.
            % the terms in g are weighted by the quaderature weight.
            sumfact{1}(k,l) = g(1) * w;
            sumfact{2}(k,l) = g(2) * w;
            sumfact{3}(k,l) = g(3) * w;
        end
    end
end