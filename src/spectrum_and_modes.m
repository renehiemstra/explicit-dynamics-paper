function [D, X] = spectrum_and_modes(A)
    
    % compute eigenvalues
    [X, D] = eig(A); X = real(X); D = real(D); 

    % sort eigenvalues and modes
    D = diag(D); [D,I] = sort(D); D = real(sqrt(D));
    X = X(:,I); 
end