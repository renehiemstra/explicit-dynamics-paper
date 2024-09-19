% perform the quadrature loop
function sumfact = mat_quadrature_loop_weighted_3d(m1, m2, m3, kappa, ...
        sumfact, ...
        weights, ...
        forcing, ...
        displacement_ders,...
        jacobian,...
        hessian)

    % for each quadrature point
    for m=1:m3
        for l=1:m2
            for k=1:m1
                % compute jacobian matrix of first and second derivatives
                % @ quadrature point (k,l,m)
                [A,DA] = compute_jacobian_and_hessian(jacobian, hessian, k, l, m);
    
                % compute the determinant of the jacobian matrix
                % @ quadrature point (k,l,m)
                c = det(A);
    
                % retreive quadrature weight
                w = weights{1}(k) * weights{2}(l) * weights{3}(m);
    
                % compute the gradient of the determinant of the jacobian matrix
                % @ quadrature point (k,l,m)
                for i=3:-1:1
                    g_grad(i) = g_grad_i(A, DA, i);
                end

                % compute the gradient of the determinant of the jacobian matrix
                % @ quadrature point (k,l,m)   
                c_grad = A' \ g_grad';

                % compute the gradient of the displacement field
                % @ quadrature point (k,l,m)
                d_grad = A' \ [displacement_ders{1}(k,l,m); displacement_ders{2}(k,l,m); displacement_ders{3}(k,l,m)];
                
                % heat flux @ quadrature point (k,l,m)
                q = kappa * d_grad;
                
                % first two elements of 'g' relate to the regular gradient.
                % the third element of 'g' relates to taking the gradient of
                % the (1/c) term and the addition of the forcing term.
                g = [-A \ q; 
                    (dot(q, c_grad) / c) + forcing(k,l,m)];
                
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

function h = g_grad_i(A, DA, i)
    h = 0;
    for j=1:3
        h = h + DA{j}(i,i) * minor(A,j,i);
        i1 = mod(i,3)+1; i2 = mod(i+1,3)+1;
        j1 = mod(j,3)+1; j2 = mod(j+1,3)+1;
        h = h + A(j,i) * (DA{j1}(i,i1) * A(j2,i2) + A(j1,i1) * DA{j2}(i,i2) -DA{j1}(i,i2) * A(j2,i1) - A(j1,i2) * DA{j2}(i,i1));
    end
end

function M = minor(A, i, j)
    M = det([A(mod(i,  3)+1, mod(j,3)+1)   A(mod(i,  3)+1, mod(j+1,3)+1)  ;
             A(mod(i+1,3)+1, mod(j,3)+1)   A(mod(i+1,3)+1, mod(j+1,3)+1)]);
end

function [A, DA] = compute_jacobian_and_hessian(dX, ddX, k, l, m)
    % get jacobian
    A = [dX{1,1}(k, l, m) dX{1,2}(k, l, m) dX{1,3}(k, l, m);
         dX{2,1}(k, l, m) dX{2,2}(k, l, m) dX{2,3}(k, l, m);
         dX{3,1}(k, l, m) dX{3,2}(k, l, m) dX{3,3}(k, l, m)];

    % get hessian in each coordinate direction
    for i=3:-1:1
        DA{i} = [ddX{1,i}(k, l, m) ddX{6,i}(k, l, m) ddX{5,i}(k, l, m);
                 ddX{6,i}(k, l, m) ddX{2,i}(k, l, m) ddX{4,i}(k, l, m);
                 ddX{5,i}(k, l, m) ddX{4,i}(k, l, m) ddX{3,i}(k, l, m)];
    end
end