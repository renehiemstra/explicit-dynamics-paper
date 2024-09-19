% perform the quadrature loop
function sumfact = mat_quadrature_loop_standard_3d(m1, m2, m3, ...
        kappa,...
        sumfact, ...
        weights, ...
        forcing, ...
        displacement_ders,...
        jacobian)

    % for each quadrature point
    for m=1:m3
        for l=1:m2
            for k=1:m1
                % compute jacobian matrix of first and second derivatives
                % @ quadrature point (k,l,m)
                A = compute_jacobian(jacobian, k, l, m);
    
                % compute the determinant of the jacobian matrix
                % @ quadrature point (k,l,m)
                c = det(A);
    
                % retreive quadrature weight
                w = weights{1}(k) * weights{2}(l) * weights{3}(m) * c;
    
                % compute the gradient of the displacement field
                % @ quadrature point (k,l)
                d_grad = A' \ [displacement_ders{1}(k,l,m); displacement_ders{2}(k,l,m); displacement_ders{3}(k,l,m)];
                
                % heat flux @ quadrature point (k,l)
                q = kappa * d_grad;
                
                % first two elements of 'g' relate to the regular gradient.
                % the third element of 'g' relates to taking the gradient of
                % the (1/c) term and the addition of the forcing term.
                g = [-A \ q; 
                    forcing(k,l,m)];
                
                % add results to sumfactorization objects.
                % the terms in g are weighted by the quaderature weight.
                sumfact{1}(k,l,m) = g(1) * w;
                sumfact{2}(k,l,m) = g(2) * w;
                sumfact{3}(k,l,m) = g(3) * w;
                sumfact{4}(k,l,m) = g(4) * w;
            end
        end
    end
end

function A = compute_jacobian(dX, k, l, m)
    A = [dX{1,1}(k, l, m) dX{1,2}(k, l, m) dX{1,3}(k, l, m);
         dX{2,1}(k, l, m) dX{2,2}(k, l, m) dX{2,3}(k, l, m);
         dX{3,1}(k, l, m) dX{3,2}(k, l, m) dX{3,3}(k, l, m)];
end