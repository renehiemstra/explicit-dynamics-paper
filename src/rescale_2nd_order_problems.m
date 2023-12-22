function X = rescale_2nd_order_problems(X, f, bcs)
    % initialize
    n = size(X,2);
    
    % recale eigenmodes
    for k=1:n
        if (strcmp(bcs,'fixed'))
            s = sign(X(2,k) - X(1,k)); 
        elseif (strcmp(bcs,'free'))
            s = sign(X(1,k));
        elseif (strcmp(bcs,'fixed-free'))
            s = sign(X(1,k));
        end
        
        X(:,k) = X(:,k) * f(X, k) * s;
    end
end

