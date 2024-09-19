function f = evaluate_field_2d(field, u, v)

    % compute basis functions
    bu = field.basis{1}.eval.trialfuns(u');
    bv = field.basis{2}.eval.trialfuns(v');

    % compute coordinates
    f = tensorconstract(field.coeffs, bu{1}, bv{1}, 2);
end