function M = system_matrix(basis, rule, i, j)
    M = basis.testfuns{i}' * (rule.weights .* basis.trialfuns{j});
end