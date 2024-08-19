function M = mass_matrix(basis, rule)
    M = basis.testfuns{1}' * (rule.weights .* basis.trialfuns{1});
end