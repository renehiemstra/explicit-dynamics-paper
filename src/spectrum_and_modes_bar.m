function [Romega, Rmodes, dtcrit]  = spectrum_and_modes_bar(method, p, ndofs, problem, remove_outliers)

    % compute spectrum and modes
    if strcmp(method, 'standard')
        [S, X, kts]  = second_order_problems_spectrum_and_modes_1d(p, ndofs, problem, remove_outliers);
    elseif strcmp(method, 'fast')
        [S, X, kts]  = second_order_problems_spectrum_and_modes_1d_fast(p, ndofs, problem, remove_outliers);
    elseif strcmp(method, 'lumped')
        [S, X, kts]  = second_order_problems_spectrum_and_modes_1d_lumped(p, ndofs, problem, remove_outliers);
    end

    % boundary conditions and analytical solution of vibrating bar
    if strcmp('free', problem)
        kappa = 1;
        f = @(x,k) cos(k.*pi.*x);
        omega = @(k) k.*pi;
    elseif strcmp('fixed', problem)
        kappa = 0;
        f = @(x,k) sin(k.*pi.*x);
        omega = @(k) k.*pi;
    elseif strcmp('fixed-free', problem)
        kappa = 0;
        f = @(x,k) cos((k-1/2) .* pi .* x);
        omega = @(k) (k-1/2) * pi;
    else
        error('Not correct definition of boundary conditions');
    end
    
    %% Compute quadrature and basisfunctions to expand the modes
    
    % compute L^2 error in the mode shapes
    [x,w] = lgwtglobal(p+1,unique(kts)); 

    % compute basis functions at quadrature points
    N = spcol(kts,p+1,x);
    W = spdiags(w, [0], length(w), length(w));
    
    %% compute error in modeshape and eigenvalue    
    m = length(S);

    % compute analytical eigenvalues
    j1 = (1+kappa:m)';

    % compute L^2 error in the mode shapes
    e = @(X,z) f(x,z-kappa) - N * X(:,z);
    l2e = @(X,z) sqrt(sum(W * e(X,z).^2, 1)) ./ sqrt(sum(W * f(x,z).^2, 1));
    
    % relative value of eigenvalue
    Romega =  S(j1) ./ omega(j1-kappa);
%     Romega =  (abs(S(j1)-omega(j1-kappa))) ./ omega(j1-kappa);
    
    % error in the mode shape
    Rmodes  = l2e(X,j1');    
    
    % compute critical timestep
    wh = max(S);
    dtcrit = 2 / wh;