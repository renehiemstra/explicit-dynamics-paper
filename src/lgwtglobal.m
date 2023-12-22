function [X,W] = lgwtglobal(n,u)

    % initialize
    m = length(u)-1;

    % Gauss-legendre points on interval [0,1]
    [x, w] = lgwt(n,0.0,1.0);
    
    % global quadrature rule
    X = zeros(n,m); W = zeros(n,m); 
    for k=1:m
        h = u(k+1)-u(k);
        X(:,k) = u(k) + h * x;
        W(:,k) = h * w;
    end
    X = X(:); W = W(:);
end

