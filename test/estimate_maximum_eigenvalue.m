function rc = estimate_maximum_eigenvalue(A, x0, maxit)

    % random input vector
    x.next = normalize(x0);
    
    % perform power iteration
    rc_current = -1.00e+15;
    it = 0;
    while true

        % normalize vector and apply operator
        x.current = normalize(x.next);
        x.next = A(x.current);
        
        % compute rayleigh coefficient
        rc = dot(x.current, x.next) / dot(x.current, x.current);
        
        % break  if errortol is achieved
        if relative_distance(rc_current, rc)< 1e-8
            break
        end
        % break if maximum number of iterations is achieved
        if it > maxit
            warning("Reached maximum number of iterations.")
            break
        end

        % initialize next iteration
        rc_current = rc;
        it = it + 1;
    end
end

% normalize input vector
function x = normalize(v)
    x = v / norm(v);
end

% compute relative distance between two numbers
function d = relative_distance(a, b)
    d = abs(b-a) / abs(b);
end