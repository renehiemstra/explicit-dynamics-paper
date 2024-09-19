function f = evaluate_field_3d(field, u, v, w)

    % compute basis functions
    bu = field.basis{1}.eval.trialfuns(u');
    bv = field.basis{2}.eval.trialfuns(v');
    bw = field.basis{3}.eval.trialfuns(w');

    % compute coordinates
    f = tensorconstract(field.coeffs, bu{1}, bv{1}, bw{1}, 2);
end