function A = constraint_matrix_periodicity(p, kts)
    A = spcol(kts, p+1, kts(1)*ones(p,1)) - spcol(kts, p+1, kts(end)*ones(p,1));
end


