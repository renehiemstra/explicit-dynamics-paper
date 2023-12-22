function A = constraint_matrix_annulus(p, kts, remove_outliers)
    A_a = apply_constraint(p,kts,kts(1),remove_outliers);
    A_b = apply_constraint(p,kts,kts(end),remove_outliers);
    A = [A_a; A_b];
end

function A = apply_constraint(p,kts,r,remove_outliers)
    d = spcol(kts, p+1, r*ones(p+1,1));
    A = d(1,:);
    if remove_outliers
        A1 = d(3,:) + (1/r) * d(2,:);
        A = [A; A1];
        if p>4
            A2 = d(5,:) + (2/r) * d(4,:) - (1/r^2) * d(3,:) + (1/r^3) * d(2,:);
            A = [A; A2];
        end
    end
end


