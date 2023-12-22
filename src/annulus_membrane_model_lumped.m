function [kts_r, kts_h, C, L, A, K] = annulus_membrane_model_lumped(p, ne_r, ne_h, remove_outliers)

    % radius of 2nd and 4th zero
    r_a = 1;
    r_b = 17.61596604980483 / 11.064709488501170;
    
    % open knot vectors
    kts_r  = [r_a*ones(1,p) linspace(r_a,r_b,ne_r+1) r_b*ones(1,p)];
    kts_h  = [0.0*ones(1,p) linspace(0.0,2*pi,ne_h+1) 2*pi*ones(1,p)];
    
    % compute extraction operator
    [E_r, E_h] = extraction_operator_annulus(p, kts_r, kts_h, remove_outliers);
    
    % compute mass and stiffness matrices
    [B_r, K_r, B_h, K_h, A_r, A_h, L_r, L_h] = univariate_mass_and_stifness_matrices(p, kts_r, kts_h, E_r, E_h);
    
    % global stiffness matrix
    K = kron(K_r,B_h) + kron(B_r,K_h); 
    A = kron(A_r, A_h);
    L = kron(L_r, L_h);
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

function [B_r, K_r, B_h, K_h, A_r, A_h, L_r, L_h] = univariate_mass_and_stifness_matrices(p, kts_r, kts_h, E_r, E_h)

    % compute L^2 error in the mode shapes
    [x_r,w_r] = lgwtglobal(2*p+1,unique(kts_r));
    [x_h,w_h] = lgwtglobal(2*p+1,unique(kts_h(p+1:end-p)));
    xx_r = repmat(x_r',2,1); xx_r = xx_r(:);
    xx_h = repmat(x_h',2,1); xx_h = xx_h(:);
    
    % functions at the quadrature points
    tmp = sparse(spcol(kts_r,p+1,xx_r));
    N_r = {tmp(1:2:end, :), tmp(2:2:end, :)};
    tmp = sparse(spcol(kts_h,p+1,xx_h));
    N_h = {tmp(1:2:end, :), tmp(2:2:end, :)};
    
    % compute mass matrices without extraction
    M_r = N_r{1}' * diag(w_r .* x_r) * N_r{1};
    M_h = N_h{1}' * diag(w_h) * N_h{1};
    
    % compute system matrices with extraction
    B_h = E_h' * N_h{1}' * diag(w_h) * N_h{1} * E_h;
    K_r = E_r' * N_r{2}' * diag(w_r .* x_r) * N_r{2} * E_r;
    B_r = E_r' * N_r{1}' * diag(w_r ./ x_r) * N_r{1} * E_r;
    K_h = E_h' * N_h{2}' * diag(w_h) * N_h{2} * E_h;
    L_r = E_r' * M_r * E_r;
    L_h = E_h' * M_h * E_h;

    % lumped inverse mass
    A_r = inv(E_r' * diag(sum(M_r,2)) * E_r);
    A_h = inv(E_h' * diag(sum(M_h,2)) * E_h);
end