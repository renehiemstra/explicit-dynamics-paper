function [S, X, kts]  = second_order_problems_spectrum_and_modes_1d_fast(p, ndofs, bcs, remove_outliers)

    % set maximum number of constraints
    q=1;
    if remove_outliers
        q = p+1;
    end
    
    % boundary conditions and analytical solution of vibrating bar
    if strcmp('free', bcs)
        v1 = 2:2:q;
        v2 = 2:2:q;
        f = @(x,k) cos(pi*x*k);
    elseif strcmp('fixed', bcs)
        v1 = 1:2:q;
        v2 = 1:2:q;
        f = @(x,k) sin(pi*x*k);
    elseif strcmp('fixed-free', bcs)
        v1 = 2:2:q;
        v2 = 1:2:q;
        f = @(x,k) cos((k-1/2)*pi*x);
    else
        error('Not correct definition of boundary conditions');
    end
    
    % initialize
    ne = ndofs - p + length(v1)+length(v2);
    kts = [zeros(1,p) linspace(0.0,1.0,ne+1) ones(1,p)]';
   
    % case with or without outliers  
    C = extraction_operator(constraint_matrix(p, kts, v1, v2));
    
    %% Compute matrices
    
    % compute L^2 error in the mode shapes
    [x,w] = lgwtglobal(p+1,unique(kts)); 
    xx = repmat(x',2,1); xx = xx(:);

    % compute basis functions at quadrature points
    N = spcol(kts,p+1,xx);
    W = spdiags(w, [0], length(w), length(w));

    % compute stiffness and right hand side
    K = N(2:2:end, :)' * W * N(2:2:end, :); KK = C' * K * C;
    M = N(1:2:end, :)' * W * N(1:2:end, :);

    % approximate inverse
    MM = inv(C' * sparse(inv(approximate_l2_inverse(p,kts,false,false))) * C);
    
    %% Compute modes and spectrum
    
    % compute normalized frequencies
    [S,X] = spectrum_and_modes(MM * KK);
    
    % rescaling function
    g = @(X,z) sqrt(f(x,z)' * W * f(x,z))   / sqrt(X(:,z)' * M * X(:,z));

    % rescale modes
    X = rescale_2nd_order_problems(C * X, g, bcs);