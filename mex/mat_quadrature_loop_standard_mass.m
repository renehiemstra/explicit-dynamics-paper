% perform the quadrature loop
function sumfact = mat_quadrature_loop_standard_mass(m1, m2, ...
        sumfact, ...
        weights, ...
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
            
            % add results to sumfactorization objects.
            % the terms in g are weighted by the quaderature weight.
            sumfact{1}(k,l) = w;
        end
    end
end