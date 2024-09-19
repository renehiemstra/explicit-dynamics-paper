% perform the quadrature loop
function sumfact = mat_quadrature_loop_standard_mass_3d(m1, m2, m3, ...
        sumfact, ...
        weights, ...
        jacobian)

    % for each quadrature point
    for m=1:m3
        for l=1:m2
            for k=1:m1
                % compute jacobian matrix of first and second derivatives
                % @ quadrature point (k,l)
                A = compute_jacobian(jacobian, k, l, m);
    
                % compute the determinant of the jacobian matrix
                % @ quadrature point (k,l)
                c = det(A);
    
                % retreive quadrature weight
                w = weights{1}(k) * weights{2}(l) * weights{3}(m) * c;
                
                % add results to sumfactorization objects.
                % the terms in g are weighted by the quadrature weight.
                sumfact{1}(k,l,m) = w;
            end
        end
    end
end

function A = compute_jacobian(dX, k, l, m)
    A = [dX{1,1}(k, l, m) dX{1,2}(k, l, m) dX{1,3}(k, l, m);
         dX{2,1}(k, l, m) dX{2,2}(k, l, m) dX{2,3}(k, l, m);
         dX{3,1}(k, l, m) dX{3,2}(k, l, m) dX{3,3}(k, l, m)];
end