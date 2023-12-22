function [kts_r, kts_h, C, L, A, K] = annulus_membrane_model_fast(p, ne_r, ne_h, remove_outliers)

    % radius of 2nd and 4th zero
    r_a = 1;
    r_b = 17.61596604980483 / 11.064709488501170;
    
    % open knot vectors
    kts_r  = [r_a*ones(1,p) linspace(r_a,r_b,ne_r+1) r_b*ones(1,p)];
    kts_h  = [0.0*ones(1,p) linspace(0.0,2*pi,ne_h+1) 2*pi*ones(1,p)];
    
    % compute extraction operator
    [E_r, E_h] = extraction_operator_annulus(p, kts_r, kts_h, remove_outliers);
    
    % compute mass and stiffness matrices
    [M_r, M_h, K_r, K_h, B_r, C_r, L_r, L_h] = univariate_mass_and_stifness_matrices(p, kts_r, kts_h, E_r, E_h);

    % global stiffness matrix
    L = kron(L_r, L_h);
    A = kron(sparse(inv(M_r)), sparse(inv(M_h)));
    K = -kron(C_r, M_h) + kron(K_r, M_h) + kron(B_r, K_h);
    C = kron(E_r, E_h);
end

function [E_r, E_h] = extraction_operator_annulus(p, kts_r, kts_h, remove_outliers)

    % Extraction operator radial direction
    A = constraint_matrix_annulus(p, kts_r, remove_outliers);
    E_r = sparse(extraction_operator(A));
    
    % extraction operator angular direction
    A = constraint_matrix_periodicity(p, kts_h);
    E_h = sparse(extraction_operator(A));
end

function [M_r, M_h, K_r, K_h, B_r, C_r, L_r, L_h] = univariate_mass_and_stifness_matrices(p, kts_r, kts_h, E_r, E_h)

    % compute L^2 error in the mode shapes
    [x_r,w_r] = lgwtglobal(2*p+1,unique(kts_r));
    [x_h,w_h] = lgwtglobal(2*p+1,unique(kts_h(p+1:end-p)));
    xx_r = repmat(x_r',2,1); xx_r = xx_r(:);
    xx_h = repmat(x_h',2,1); xx_h = xx_h(:);
    
    % compute basis functions at quadrature points
    N_r = spcol(kts_r,p+1,xx_r) * E_r;
    N_h = spcol(kts_h,p+1,xx_h) * E_h;
    
    % compute mass matrices
%     M_r = sparse(N_r(1:2:end, :)' * diag(w_r) * N_r(1:2:end, :));
%     M_h = sparse(N_h(1:2:end, :)' * diag(w_h) * N_h(1:2:end, :));
    M_r = E_r' * sparse(inv(approximate_l2_inverse(p,kts_r,false,false))) * E_r;
    M_h = E_h' * sparse(inv(approximate_l2_inverse(p,kts_h,false,false))) * E_h;
    K_r = sparse(N_r(2:2:end, :)' * diag(w_r) * N_r(2:2:end, :));
    K_h = sparse(N_h(2:2:end, :)' * diag(w_h) * N_h(2:2:end, :));
    B_r = sparse(N_r(1:2:end, :)' * diag(w_r ./ x_r.^2) * N_r(1:2:end, :));
    C_r = sparse(N_r(1:2:end, :)' * diag(w_r ./ x_r) * N_r(2:2:end, :));
    L_r = sparse(N_r(1:2:end, :)' * diag(w_r .* x_r) * N_r(1:2:end, :));
    L_h = sparse(N_h(1:2:end, :)' * diag(w_h) * N_h(1:2:end, :));
end